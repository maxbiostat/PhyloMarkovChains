% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighbourhood.r
\name{neighbourhood_size}
\alias{neighbourhood_size}
\title{Compute the size of rSPR neighbourhood}
\usage{
neighbourhood_size(phy, method = "song")
}
\arguments{
\item{phy}{Phylogenetic Tree}

\item{method}{Method used, could be "song" or "matsen"
where \code{song} uses the method by Song (2003) and \code{matsen} uses the
approach in Whidden and Matsen IV (2017). Default is method = "song".}
}
\value{
The size of the neighbourhood of a phylogenetic tree in the rSPR Graph
}
\description{
Compute the size of rSPR neighbourhood
}
\examples{
t <- ape::rtree(5, rooted = TRUE)
neighbourhood_size(t, method = "song")
neighbourhood_size(t, method = "matsen")


}
\references{
Song, Y. S. (2003). On the combinatorics of rooted binary
phylogenetic trees. \emph{ Annals of Com-binatorics, 7(3):365?37}

Whidden, C. and Matsen IV, F. A. (2017). Ricci-Ollivier
curvature of the rooted phylogenetic subtree-prune-regraft graph.
\emph{ Theoretical Computer Science, 699:1?2}
}
