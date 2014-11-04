setClass(Class = "coalescent",
	representation(	node = "numeric",
					ancestor = "numeric",
					label = "character",
					branch.length = "numeric")
)

setClass(Class = "history",
	representation(state = "numeric"),
	contains = "coalescent"
)

is.coalescent <- function (x) 
{
    res <- (is(x,"coalescent") & validObject(x))
    return(res)
}

is.history <- function (x) 
{
    res <- (is(x,"history") & validObject(x))
    return(res)
}

as.history <- function(tree) {
	if(!is.coalescent(tree)) {
		stop(paste(tree," is not a valid coalescent object...",sep = ""))
	}	
	tree.with.mutations <- new("history")
	tree.with.mutations@node <- tree@node
	tree.with.mutations@ancestor <- tree@ancestor
	tree.with.mutations@label <- tree@label
	tree.with.mutations@branch.length <- tree@branch.length
	return(tree.with.mutations)	
}

SimulateCoalescentTree <- function(method,sample,current,ancestral,time) {
	.C('SimulateCoalescentTree',
	as.integer(method),
	as.integer(sample),
	as.integer(current),
	as.integer(ancestral),
	as.integer(time),
	PACKAGE = 'coalesceR')
}

read.tree <- function(infile = "",tree) {											# parse a tree in Newick format
	if (infile != "") {
		if(!file.exists(infile)) {
			stop(paste("The file ",infile," does not exist. Check the input file name...",sep = ""))
		}
		tree <- scan(infile,sep = "\n",what = "character",quiet = TRUE)
	}
	newick.tree <- unlist(strsplit(tree,split = ""))								# Split the string tree into single characters
	if (newick.tree[length(newick.tree)] != ";") {
		stop("Problem with Newick format: it must end with a ';'")
	}
	if (length(grep("\\(",newick.tree)) != length(grep(")",newick.tree))) {
		stop("Problem with Newick format: the numbers of '(' and ')' differ")
	}
	sep <- sort(c(grep("\\(",newick.tree),grep(")",newick.tree),grep(":",newick.tree),grep(",",newick.tree), grep(";",newick.tree)))
	tmp <- {}																		# group the labels and branch lengths into a single item
	for (i in 1:(length(sep) - 1)) {
		if (length(sep[i]:sep[i + 1]) > 2) {
			tmp <- append(tmp,newick.tree[sep[i]])
			tmp <- append(tmp,paste(newick.tree[(sep[i] + 1):(sep[i + 1] - 1)],collapse = ""))
		} else {
			tmp <- append(tmp,newick.tree[sep[i]])
		}
	}	
	newick.tree <- append(tmp,newick.tree[length(newick.tree)])
	leafs <- length(grep(",",newick.tree)) + 1
	nodes <- length(grep("\\(",newick.tree))
	tree <- new("coalescent")
	tree@node <- seq(1,(nodes + leafs - 1))
	i <- 1																			# Need to parse successive items
	current.node <- (nodes + leafs)													# Start at the MRCA
	current.leaf <- 1 																# Count the tips...
	next.node <- current.node - 1													# 
	while (newick.tree[i] != ";") {													# 
		if (newick.tree[i] == "(") {												# Create a daughter node
			i <- i + 1
			if (is.na(match(newick.tree[i],c(",","(",")",":",";")))) {				# The next node is a leaf
				tree@ancestor[current.leaf] <- current.node							# Record the ancestor of the leaf
				tree@label[current.leaf] <- newick.tree[i]							# Record the label
				i <- i + 1
				if (newick.tree[i] == ":") {										# Is there a branch length?
					i <- i + 1
					tree@branch.length[current.leaf] <- as.numeric(newick.tree[i])	# Record the branch length
					i <- i + 1
				}
				current.leaf <- current.leaf + 1									# Increment the leaf
			} else if (newick.tree[i] == "(") {										# The next node is an internal node
				tree@ancestor[next.node] <- current.node
				current.node <- next.node
				next.node <- next.node - 1
			}
		} else if (newick.tree[i] == ")") {
			i <- i + 1
			if (newick.tree[i] == ":") {											# Is there a branch length?
				i <- i + 1
				tree@branch.length[current.node] <- as.numeric(newick.tree[i])		# Record the branch length
				i <- i + 1
			}
			if (current.node < (nodes + leafs)) {									# ???!!!!!
				current.node <- tree@ancestor[current.node]							# Go back to the ancestral node
			}
		} else if (newick.tree[i] == ",") {
			i <- i + 1
			if (is.na(match(newick.tree[i],c(",","(",")",":",";")))) {				# The next node is a leaf
				tree@ancestor[current.leaf] <- current.node							# Record the ancestor of the leaf
				tree@label[current.leaf] <- newick.tree[i]							# Record the label
				i <- i + 1
				if (newick.tree[i] == ":") {										# Is there a branch length?
					i <- i + 1
					tree@branch.length[current.leaf] <- as.numeric(newick.tree[i])	# Record the branch length
					i <- i + 1
				}
				current.leaf <- current.leaf + 1									# Increment the leaf
			} else if (newick.tree[i] == "(") {										# The next node is an internal node
				tree@ancestor[next.node] <- current.node
				current.node <- next.node
				next.node <- next.node - 1
			}
		}
	}
	return(tree)
}

sim.tree <- function(method,sample,current,ancestral,time) {
	if (method == "generations") {
		SimulateCoalescentTree(method = 0,sample = sample,current = current,ancestral = ancestral,time = time)
	} else if (method == "hudson") {
		SimulateCoalescentTree(method = 1,sample = sample,current = current,ancestral = ancestral,time = time)
	}
	tree <- read.tree(infile = "tree.dat")
	file.remove("tree.dat")
	return(tree)
}

draw.tree <- function(tree,labels = TRUE){	
	if(!is.coalescent(tree)) {
		stop(paste(tree," is not a valid coalescent object...",sep = ""))
	}	
	internal.nodes <- sort(unique(tree@ancestor))									# Get the list of internal nodes	
	leafs <- tree@node[1:length(tree@label)]										# Get the list of leafs
	nodes <- tree@node[-leafs]														# Get the list of nodes (but the mrca)
	mrca <- max(internal.nodes)
	x <- vector(mode = "numeric",length = (length(leafs) + length(internal.nodes)))
	x[1:length(leafs)] <- leafs														# X-coordinates of the leafs, by definition
	for (i in internal.nodes) {
		x[i] <- mean(x[tree@node[tree@ancestor == i]])								# Compute the x-coordinates of the internal nodes
	}
	y <- vector(mode = "numeric",length = (length(leafs) + length(internal.nodes)))	# Compute the y-coordinates of all the nodes
	for (i in c(leafs,nodes)) {														# For each lineage, go back to the mrca
		idx <- i
		while (idx != mrca) {
			y[i] <- y[i] + tree@branch.length[idx]									# Compute the height of each lineage
			idx <- tree@ancestor[idx]
		}	
	}
	y <- max(y) - y
	x.scaling <- (max(x[leafs]) - min(x[leafs])) * 0.1
	y.scaling <- (max(y) - min(y)) * 0.1
	dev.new()
	plot(x[leafs],y[leafs],xlim = c(min(x[leafs]) - x.scaling,max(x[leafs]) + x.scaling),ylim = c(min(y) - y.scaling,max(y) + y.scaling),frame.plot = FALSE,xaxt = "n",xlab = "",ylab = "Coalescence time",pch = 16)
	if (labels) {
		text(x[leafs],y[leafs], tree@label[leafs],pos = 1)							# Plot the labels of all the leafs
	}
	for (i in sort(internal.nodes,decreasing = TRUE)) {								# Plot the coalescent tree
		lineages <- tree@node[tree@ancestor == i]									# List of the lineages that descend from i
		for (j in lineages) {
			lines(c(x[j],x[j]),c(y[j],y[i]))										# Plot the vertical lines for each descendant
		}
		lines(c(min(x[lineages]),max(x[lineages])),c(y[i],y[i]))					# Plot the horizontal line for the ancestral lineage
	}
	points(x[mrca],y[mrca],pch = "|")
	text(x[mrca],y[mrca],"MRCA",pos = 3)
}

add.mutations <- function(tree,mu,graphics = TRUE,labels = TRUE) {
	if(!is.coalescent(tree)) {
		stop(paste(tree," is not a valid coalescent object...",sep = ""))
	}	
	tree.with.mutations <- as.history(tree)
	leafs <- tree@node[1:length(tree@label)]										# Get the list of leafs
	nodes <- tree@node[-leafs]														# Get the list of nodes (but the mrca)
	internal.nodes <- sort(unique(tree@ancestor))
	mrca <- max(internal.nodes)
	if (graphics) {
		x <- vector(mode = "numeric",length = (length(leafs) + length(internal.nodes)))
		x[1:length(leafs)] <- leafs													# X-coordinates of the leafs, by definition
		for (i in internal.nodes) {
			x[i] <- mean(x[tree@node[tree@ancestor == i]])							# Compute the x-coordinates of the internal nodes
		}
		y <- vector(mode = "numeric",length = (length(leafs) + length(internal.nodes)))	# Compute the y-coordinates of all the nodes
		for (i in c(leafs,nodes)) {													# For each lineage, go back to the mrca
			idx <- i
			while (idx != mrca) {
				y[i] <- y[i] + tree@branch.length[idx]								# Compute the height of each lineage
				idx <- as.numeric(tree@ancestor[idx])
			}	
		}
		y <- max(y) - y		
		x.scaling <- (max(x[leafs]) - min(x[leafs])) * 0.1
		y.scaling <- (max(y) - min(y)) * 0.1
		dev.new()
		plot(x[leafs],y[leafs],xlim = c(min(x[leafs]) - x.scaling,max(x[leafs]) + x.scaling),ylim = c(min(y) - y.scaling,max(y) + y.scaling),frame.plot = FALSE,xaxt = "n",xlab = "",ylab = "Coalescence time",pch = 16)
		if (labels) {
			text(x[leafs],y[leafs], tree@label[leafs],pos = 1)						# Plot the labels of all the leafs
		}
		for (i in sort(internal.nodes,decreasing = TRUE)) {							# Plot the coalescent tree
			lineages <- tree@node[tree@ancestor == i]								# List of the lineages that descend from i
			for (j in lineages) {
				lines(c(x[j],x[j]),c(y[j],y[i]))									# Plot the vertical lines for each descendant
			}
			lines(c(min(x[lineages]),max(x[lineages])),c(y[i],y[i]))				# Plot the horizontal line for the ancestral lineage
		}
		points(x[mrca],y[mrca],pch = "|")
		text(x[mrca],y[mrca],"MRCA",pos = 3)
	}
	iam <- 0
	for (i in sort(internal.nodes,decreasing = TRUE)) {								# Plot the coalescent tree
		lineages <- tree@node[tree@ancestor == i]									# List of the lineages that descend from i
		for (j in lineages) {
			nbr.mutations <- rpois(1,(tree@branch.length[j] * mu))
			if (nbr.mutations > 0) {
				for (k in 1:nbr.mutations) {
					iam <- iam + 1
				}
				tree.with.mutations@state[j] <- iam
			} else {
				if (tree@ancestor[j] < mrca) {
					tree.with.mutations@state[j] <- tree.with.mutations@state[i]
				} else {
					tree.with.mutations@state[j] <- 0
				}				
			}
			if (graphics) {		
				y.mut <- y[i] - runif(nbr.mutations) * tree@branch.length[j]
				x.mut <- rep(x[j], nbr.mutations)
				points(x.mut,y.mut,col = "red",pch = "-")
			}
		}
	}
	
	if (graphics) {		
		sample <- tree.with.mutations@state[leafs]
		alleles <- unique (sample)
		color.alleles <- rainbow(length(alleles))
		names(color.alleles) <- alleles
		points(x[leafs],y[leafs],pch = 16,col = color.alleles[as.character(sample)])
	}	
	return(tree.with.mutations)
}

heterozygosity <- function(tree.with.mutations){	
	if(!is.history(tree.with.mutations)) {
		stop(paste(tree.with.mutations," is not a valid history object...",sep = ""))
	}	
	leafs <- tree.with.mutations@node[1:length(tree.with.mutations@label)]			# Get the list of leafs
	sample <- table(tree.with.mutations@state[leafs])
	n <- sum(sample)
	frequency <- as.matrix(sample) / n
	H <- (1 - sum(frequency^2)) * (n / (n -1))
	return(H)
}

number.of.alleles <- function(tree.with.mutations){
	if(!is.history(tree.with.mutations)) {
		stop(paste(tree.with.mutations," is not a valid history object...",sep = ""))
	}	
	leafs <- tree.with.mutations@node[1:length(tree.with.mutations@label)]			# Get the list of leafs
	sample <- tree.with.mutations@state[leafs]
	k <- unique(sample)
	return(length(k))	
}

TMRCA <- function(tree) {
	if(!is.coalescent(tree)) {
		stop(paste(tree," is not a valid coalescent or history object...",sep = ""))
	}
	internal.nodes <- sort(unique(tree@ancestor))									# Get the list of internal nodes	
	mrca <- max(internal.nodes)
	tmrca <- 0
	j = 10
	repeat {
		tmrca <- tmrca + tree@branch.length[j]
		j <- tree@ancestor[j]
		if (j >= mrca) break
	}
	return(tmrca)
} 






