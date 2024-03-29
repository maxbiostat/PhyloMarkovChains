% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_t_mix.r
\name{t_mix_bounds}
\alias{t_mix_bounds}
\title{Title Bounds_for_t_mix_of_a_MC}
\usage{
t_mix_bounds(Tr_mat, varepislon = 0.01, pi_min)
}
\arguments{
\item{Tr_mat}{The transition matrix of the Markov Chain.}

\item{varepislon}{The maximum distance fo the Markov Chain at a given time
its stationary distribution.}

\item{pi_min}{The minimum value in the stationary distribution.}
}
\value{
An array such that in the first position a lower bound and the
second an upper bound for the mixing time of the Markov Chain with matrix
transition \code{Tr_mat}.
}
\description{
This function gives a lower and upper bounds for the mixing time of
a Markov Chain.
}
\examples{
TM <- matrix(c(0, 1/3, 1/4, 0, 0, 1/2, 0, 1/4, 1/2, 0, 1/2, 1/3, 0, 1/2, 1, 
0, 1/3, 1/4, 0, 0, 0, 0, 1/4, 0, 0), nrow = 5, ncol = 5)
p_m <- 1/12
t_mix_bounds(Tr_mat = TM, pi_min = p_m)
# Output: [1] 0.000000 7.090077

}
\references{
Levin, D. A. and Peres, Y. (2017). \emph{Markov chains and mixing times,
volume 107. American Mathematical Soc.}
}
