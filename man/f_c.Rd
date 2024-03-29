% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_c.r
\name{f_c}
\alias{f_c}
\title{Excise a clade from a tree}
\usage{
f_c(phy, clade)
}
\arguments{
\item{phy}{a phylo object.}

\item{clade}{a list of tip names in a clade.}
}
\value{
a list with two phylo objects: one is the tree without the clade and
the other is the subtree with the clade.
}
\description{
Excise a clade from a tree
}
