% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shape_lumping_aux.r
\name{get_shapes}
\alias{get_shapes}
\title{Get the shapes of all trees in a set}
\usage{
get_shapes(all_trees, ncores = 6)
}
\arguments{
\item{all_trees}{the set of trees to be compared.
Usually the output of phangorn::allTrees()}

\item{ncores}{number of cores to be used. Default is \code{ncores=6}.}
}
\value{
a binary matrix which is 1 if two trees have the same shape and 0 otherwise
}
\description{
Get the shapes of all trees in a set
}
