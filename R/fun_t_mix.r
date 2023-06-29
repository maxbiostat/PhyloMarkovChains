#' Title Bounds_for_t_mix_of_a_MC
#'
#'  This function gives a lower and upper bounds for the mixing time of 
#'  a Markov Chain.
#'  
#' @param Tr_mat The transition matrix of the Markov Chain.
#' @param varepislon The maximum distance fo the Markov Chain at a given time
#' its stationary distribution.
#' @param pi_min The minimum value in the stationary distribution.  
#'  
#' @return An array such that in the first position a lower bound and the 
#' second an upper bound for the mixing time of the Markov Chain with matrix 
#' transition \code{Tr_mat}.
#' @export t_mix_bounds
#'
#' @examples
#' TM <- matrix(c(0, 1/3, 1/4, 0, 0, 1/2, 0, 1/4, 1/2, 0, 1/2, 1/3, 0, 1/2, 1, 
#' 0, 1/3, 1/4, 0, 0, 0, 0, 1/4, 0, 0), nrow = 5, ncol = 5)
#' p_m <- 1/12
#' t_mix_bounds(Tr_mat = TM, pi_min = p_m)
#' # Output: [1] 0.000000 7.090077
#' 
#' @references 
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov chains and mixing times, 
#' volume 107. American Mathematical Soc.}
#' 
#' 
t_mix_bounds <-function(Tr_mat, varepislon = 0.01, pi_min){

  eigenva <- eigen(Tr_mat)
  max_eigenva <- max(eigenva$values[2], abs(eigenva$values[nrow(Tr_mat)]))
  t_rel <- 1/(1 - max_eigenva)
  
  upper_t_mix <- log(1/(varepislon*pi_min))*t_rel
  lower_t_mix <- (t_rel - 1)*log(1/(2*varepislon))
  
  return(c(lower_t_mix, upper_t_mix))
}





