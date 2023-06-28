
#' Get clade index
#'
#' @param x the clade as a vector of tip labels
#' @param all.splits a list of all clades (splits)
#'
#' @return an index
#' @export clade_index
clade_index <- function(x, all.splits){
  charmatch(x, all.splits)
} 

#' Construct clade indicators
#'
#' @param ind indices that should be \code{1}
#' @param L size of the indicator vector
#'
#' @return a binary vector with which clades are present (1) and which are absent (0)
#' @export clade_indicators
#'
clade_indicators <- function(ind, L){
  out <- rep(0, L)
  out[ind] <- 1
  return(out)
}
#' Get the clades in a tree
#'
#' @param tree a phylo object
#'
#' @return a vector with the 2n-1 clades
#' @export get_clades
#'
get_clades <- function(tree){
  splits <- phangorn::as.splits(tree)
  sorted <- sort(tree$tip.label)
  label.orders <- sapply(tree$tip.label, function(y) match(y, sorted))
  ordered.splits <- sapply(splits, function(x) label.orders[x])
  out <- sapply(ordered.splits,  function(x) paste(sort(x), collapse = ""))
  return(out)
} 
