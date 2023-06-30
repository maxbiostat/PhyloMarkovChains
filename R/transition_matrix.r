#' Create the Metropolis-Hastings transition matrix
#'
#' @param incidence.mat an incidence matrix (binary entries) of the SPR graph
#'
#' @return a MH transition matrix
#' @export make_MH_matrix
#'
make_MH_matrix <- function(incidence.mat) {
  K <- nrow(incidence.mat)
  if (ncol(incidence.mat) != K)
    stop("Incidence matrix is not square")
  
  neighbourhood.sizes <- colSums(incidence.mat)
  
  neigh.ratios <- matrix(NA, nrow = K, ncol = K)
  diag(neigh.ratios) <- 0
  for (i in 1:K) {
    neigh.ratios[i,] <-
      sapply(neighbourhood.sizes[i] / neighbourhood.sizes,
             function(x)
               min(1, x)) /
      neighbourhood.sizes[i]
  }
  trans_Mat <- neigh.ratios * incidence.mat
  diag(trans_Mat) <- 1 - rowSums(trans_Mat)
  
  return(trans_Mat)
}

#' 'Lazify' a Metropolis-Hastings transition matrix
#'
#' @param MHMat original MH transition matrix
#' @param rho stickyness coefficient. The probability of staying put.
#'
#' @return a MH transition matrix of a rho-lazy MH process
#' @export make_lazy_MH_matrix
#'
make_lazy_MH_matrix <- function(MHMat, rho) {
  K <- nrow(MHMat)
  if (ncol(MHMat) != K)
    stop("MH Transition matrix is not square")
  
  if (rho < 0 || rho > 1)
    stop("Rho must be a probability")
  
  trans.mat.LazyMH <- (1 - rho) * MHMat + rho * diag(1, K)
  
  return(trans.mat.LazyMH)
}
