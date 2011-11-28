#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>

#define MAXCOAL             50
#define TREEFILE           "tree.dat"

struct Node {
	int number;
	double time;                                                                      // age of the node (time = 0 at sample time)
	int allele;                                                                    // allelic state of the node
	int nbr_descendants;                                                           // number of descendants of the nodes
	int parent_nbr;                                                                // indicator variable that gives the "parent" (a number) for a node (this is for the coalescence algorithm)
	struct Node *ancestor;                                                         // address of the ancestor of the node
	struct Node *descendant[MAXCOAL];                                              // addresses of the descendants of the nodes
};

typedef struct {
	double method;
	int n;
	double time;
	double current;
	double ancestral;
	
} Parameters;

void BuildTreeGenerations(Parameters P);
void BuildTreeHudson(Parameters P);
void PrintTree(struct Node *node,FILE *treefile);
unsigned long InitializeSeed(void);

unsigned long seed;
struct Node *tree;
struct Node **list;
int nbr_nodes;

void SimulateCoalescentTree(int *method,int *sample,int *current,int *ancestral,int *time) {
	int max_nodes;
	Parameters P;  
	FILE *treefile;

	GetRNGstate();	
	P.method = *method;
	P.n = *sample;
	P.current = *current;
	P.ancestral = *ancestral;
	P.time = *time;
	max_nodes = 2 * P.n - 1;                                        // This is the maximum number of nodes in the complete genealogy
	list = (struct Node **) malloc (P.n * sizeof (struct Node *));  // These are the pointers to the remaining lineages
	tree = (struct Node *) malloc (max_nodes * sizeof(struct Node));               // This contains the full genealogy	
	treefile = fopen(TREEFILE,"w");
	if (P.method == 0) {
		BuildTreeGenerations(P);                                                            // Give a sample of genes from a (generation-by-generation) coalescent alogrithm with mutations
		PrintTree(&tree[nbr_nodes - 1],treefile);
	} else if (P.method == 1) {
		BuildTreeHudson(P);                                                            // Give a sample of genes from a (generation-by-generation) coalescent alogrithm with mutations
		PrintTree(&tree[2 * P.n - 2],treefile);
	}
	fclose(treefile);  
	free(tree);                                                                    // free the memory allocated for the structures 'tree', 'sample', and 'list'
	free(list);
	PutRNGstate();
}

void BuildTreeGenerations(Parameters P)

{
	int i,j,k;
	int sampled_nodes;
	int time = 0;
	int N;
	
	sampled_nodes = 2 * P.n - 1;                                    // Initialize the parameters for the coalescent tree
	for (i = 0; i < sampled_nodes; ++i) {                                          // Initialization of terminal nodes:
		tree[i].ancestor = NULL;                                                     // -no ancestor
		tree[i].nbr_descendants = 0;                                                 // -no descendants
		for (j = 0; j < MAXCOAL; ++j) {                                              // -no address to descendants
			tree[i].descendant[j] = NULL;
		}
	}
	for (i = 0; i < P.n; ++i) {
		tree[i].number = i;
		tree[i].time = 0.0;                                                            // 'tree[i]' [for i = 0, 1, ..., (n - 1)] are the sampled (terminal) nodes
		list[i] = tree + i;                                                          // 'list' points to the ancestral nodes of the sample
	}
	nbr_nodes = i;                                                                 // 'i' is the total number of lineages at sample time, 'nbr_nodes' is a counter for the internal nodes
	while (i > 1) {                                                                // Realization of the coalescent process for a sample of genes : the process lasts until MRCA is reached
		time++;                                                                      // Increment the number of generations
		if (time >= P.time) {
			N = P.ancestral;			
		} else {
			N = P.current;
		}
		for(j = 0; j < i; ++j) {                                                     // Attribute a "parent" to each lineage (this is to determine which lineages share the same parent, i.e. coalesce OR disperse jointly)
			list[j] -> parent_nbr = (int) (unif_rand() * N);
		}
		for(j = 0; j < (i - 1); ++j) {                                               // Loop over all pairs (j, k > j) of lineages
			if (i > 1) {                                                               // Check that the MRCA is not reached yet
				for(k = (j + 1); k < i; ++k) {
					if (list[j] -> parent_nbr == list[k] -> parent_nbr) {                  // ...then these two lineages coalesce
						if (list[j] -> time != (double) time) {                                       // The jth lineage has not coalesced yet: coalescence of a pair of lineage
							tree[nbr_nodes].parent_nbr = list[j] -> parent_nbr;                // tree[i] for i = n, (n + 1), ... are the internal nodes. The new internal node's 'parent_nbr' takes the value of list[j] -> parent_nbr
							tree[nbr_nodes].time = (double) time;                                       // The new internal node's time is the current time
							list[j] -> ancestor = tree + nbr_nodes;                            // tree[j] (to which list[j] points) now has ancestor tree + nbr_nodes (which is the address of tree[nbr_nodes])
							list[k] -> ancestor = tree + nbr_nodes;                            // tree[k] (to which list[k] points) now has ancestor tree + nbr_nodes (which is the address of tree[nbr_nodes])
							tree[nbr_nodes].nbr_descendants = 2;                               // The new internal node has two descendants (pairwise coalescence)
							tree[nbr_nodes].descendant[0] = list[j];                           // tree[nbr_nodes]'s descendant[0] is the address of list[j]
							tree[nbr_nodes].descendant[1] = list[k];                           // tree[nbr_nodes]'s descendant[1] is the address of list[k]
							tree[nbr_nodes].number = nbr_nodes;
							list[j] = tree + nbr_nodes;                                        // list[j] points to the new internal node
							list[k] = list[--i];                                               // The number (i) of lineages decreases AND THEN list[k] points to the last lineage
							++nbr_nodes;
							--k;                                                               // This is to ensure that all pairs of lineages are considered
						}
						else {                                                               // The jth lineage has coalesced already: multiple coalescence (>2 lineages)
							list[j] -> descendant[list[j] -> nbr_descendants] = list[k];       // The node to which list[j] points has one more descendant: list[k]
							list[j] -> nbr_descendants++;                                      // Increment the number of descendants for the node to which list[j] points
							list[k] -> ancestor = list[j];                                     // The node to which list[k] points has list[j] as an ancestor
							list[k] = list[--i];                                               // The number (i) of lineages decreases AND THEN list[k] points to the last lineage
							--k;                                                               // This is to ensure that all pairs of lineages are considered
						}
					}
				}
			}
		}                                                                            // This is to chose homolog lineages
	}
}

void BuildTreeHudson(Parameters P)

{
	int i,j,k;
	int sampled_nodes;
	double time = 0;
	double N;
	
	sampled_nodes = 2 * P.n - 1;                                    // Initialize the parameters for the coalescent tree
	for (i = 0; i < sampled_nodes; ++i) {                                          // Initialization of terminal nodes:
		tree[i].ancestor = NULL;                                                     // -no ancestor
		tree[i].nbr_descendants = 0;                                                 // -no descendants
		for (j = 0; j < 2; ++j) {                                              // -no address to descendants
			tree[i].descendant[j] = NULL;
		}
	}
	for (i = 0; i < P.n; ++i) {
		tree[i].number = i;
		tree[i].time = 0.0;
		list[i] = tree + i;
	}
	N = (double) P.current;
	while (i > 1) {
		time += N * (-2.0 * log(1 - unif_rand()) / (((double)i) * (i - 1)));
		if (time < (double) P.time) {
			j = (int) (unif_rand() * i);
			do {
				k = (int) (unif_rand() * i);
			}
			while (k == j);
			tree[2 * P.n - i].time = time;
			list[j] -> ancestor = tree + 2 * P.n - i;
			list[k] -> ancestor = tree + 2 * P.n - i;
			tree[2 * P.n - i].nbr_descendants = 2;
			tree[2 * P.n - i].descendant[0] = list[j];
			tree[2 * P.n - i].descendant[1] = list[k];
			tree[2 * P.n - i].number = 2 * P.n - i;
			list[j] = tree + 2 * P.n - i;
			list[k] = list[--i];
		} else {
			time = (double) P.time;
			N = (double) P.ancestral;
			P.time = 1.7976931348623157e+308;
		}
	}
}

void PrintTree(struct Node *node,
               FILE *treefile)

{
	int i;
	
	if (node -> nbr_descendants == 0) {
		fprintf(treefile,"%d:%10.6f",node -> number,(node -> ancestor) -> time);
	}
	else {
		fprintf(treefile,"(");
		for (i = 0; i < node -> nbr_descendants; ++i) {
			PrintTree(node -> descendant[i],treefile);
			if (i < ((node -> nbr_descendants) - 1)) {
				fprintf(treefile,",");
			}
		}
		if (node -> ancestor != NULL) {
			fprintf(treefile,"):%10.6f",((node -> ancestor) -> time - node -> time));
		}
		else {
			fprintf(treefile,");\n");
		}
	}
}
