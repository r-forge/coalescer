\name{simulate.tree}
\Rdversion{0.1}
\alias{simulate.tree}
\title{Simulate a Coalescent Tree}
\description{
Simulate a coalescent tree.
}
\usage{
  simulate.tree(method,sample,current,ancestral,time)
}
\arguments{
  \item{method}{simulation method: \sQuote{hudson} for the continuous time approximations; \sQuote{generations} for the generation-by-generation algorithm}
  \item{sample}{number of sampled lineages}
  \item{current}{current population size}
  \item{ancestral}{ancestral population size}
  \item{time}{time at which the population changed in size}
}
\details{
A coalescent tree is built using either the \sQuote{hudson}'s approximation in continuous time, or a generation-by-generation algorithm. A single population is simulated, and \sQuote{sample} lineages are simulated. The population is assumed to have a current effective size (\sQuote{current}), and an ancestral effective size (\sQuote{ancestral}). The change in population size occurred \sQuote{time} generations in the past
}
\value{
an object of class \sQuote{coalescent} which is a list with components:
\item{node}{the list of nodes (internal and leafs) in the tree}
\item{ancestor}{the ancestor of each node}
\item{label}{the label of each leaf}
\item{branch.length}{the length of the branch above each node}
}
\examples{
## This is to simulate a coalescent tree using the Hudson's algorithm, 
## with 20 sampled lineages. In this example, the current population size
## is 100. It was 1000 before a bottleneck that occurred 50 generations ago

tree <- simulate.tree(method = "hudson",sample = 20,current = 100,ancestral = 1000,time = 50)
}
