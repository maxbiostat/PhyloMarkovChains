clade_index <- function(x, all.splits) charmatch(x, all.splits)
clade_indicators <- function(ind, L){
  out <- rep(0, L)
  out[ind] <- 1
  return(out)
}
#' Get the clades in a tree
#'
#' @param tree a phylo object
#' @param labels labels for the tips
#'
#' @return a vector with the 2n-1 clades
#' @export get_clades
#'
get_clades <- function(tree, labels){
  splits <- phangorn::as.splits(tree, labels = labels)
  sorted <- sort(tree$tip.label)
  label.orders <- sapply(tree$tip.label, function(y) match(y, sorted))
  ordered.splits <- sapply(splits, function(x) label.orders[x])
  out <- sapply(ordered.splits,  function(x) paste(sort(x), collapse = ""))
  return(out)
} 
