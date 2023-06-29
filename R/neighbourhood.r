#' Make a ladder tree
#'
#' @param ntaxa numerical. An integer giving the number of taxa
#'
#' @return a phylo object containing a ladder treee
#' @export makeLadder
#'
makeLadder <- function(ntaxa){
  shapeTree <- function(tree){
    newTree <- tree
    nt <- ape::Ntip(newTree)
    tips <- 1:nt
    for(tip in tips){
      parent <- tree$edge[which(tree$edge[, 2] == tip), ][1]
      if(parent == ape::Ntip(tree) + 1){
        newTree$edge.length[which(newTree$edge[, 2] == tip)] <- 1
        newTree$edge.length[which(newTree$edge[, 2] == parent)] <-  newTree$edge.length[which(tree$edge[, 2] == tip)]
      }else{
        newTree$edge.length[which(newTree$edge[, 2] == tip)] <- newTree$edge.length[which(tree$edge[, 2] == parent)]
        newTree$edge.length[which(newTree$edge[, 2] == parent)] <-  newTree$edge.length[which(tree$edge[, 2] == tip)]
      }
    }
    return(newTree)
  }
  Q <- function(n, i) (i==1)
  tmp <- apTreeshape::as.phylo.treeshape(apTreeshape::rtreeshape(n = 1,
                                                    tip.number = ntaxa,
                                                    FUN = Q)[[1]])
  tmplad <- shapeTree(tmp)
  return(tmplad)
}


#' Compute the gamma of an internal node
#'
#'  This function gives the gamma of the internal node in a phylogenetic tree.
#'
#' @param phy Phylogenetic tree.
#' @param no Internal node.
#'
#' @return The corresponding gamma, defined in Song (2003), of the internal node
#' of the phylogenetic tree.
#' @export gamma_node
#'
#' @examples
#' 
#' t <- ape::rtree(5, rooted = TRUE)
#' gamma_node(t, 8)
#' # Output:  [1] 1
#'
#' @references Song, Y. S. (2003). On the combinatorics of rooted binary
#' phylogenetic trees. \emph{ Annals of Combinatorics, 7(3):365?37}
#'
#'
gamma_node <- function(phy, no) {
  A <- phangorn::Ancestors(phy, no, type = "all")
  G_no <- length(A) - 1
  return(G_no)
}

#' Sum of gammas of a phylogenetic tree
#'
#' This function gives the sum of the gammas from a phylogenetic tree.
#'
#'
#' @param phy Phylogenetic tree.
#'
#' @return The sum of the gammas from the internal nodes of a phylogenetic tree,
#' defined in Song (2003).
#' @export compute_gamma
#'
#' @examples
#' t <- ape::rtree(5, rooted = TRUE)
#' compute_gamma(t)
#' # Output:  [1] 1
#'
#' @references Song, Y. S. (2003). On the combinatorics of rooted binary
#' phylogenetic trees. \emph{ Annals of Combinatorics, 7(3):365?37}
#'
compute_gamma <- function(phy) {
  n <- length(phy$tip.label)
  N <- phy$Nnode
  q_gammas <- N - 1
  Gammas <- array(NA, dim = q_gammas)
  for (i in (n + 2):(n + N)) {
    Gammas[i - n - 1] <- gamma_node(phy, i)
  }
  return(sum(Gammas))
}

#' Size of the rSPR neighbourhood using the method of Song (2003)
#'
#' Compute the number of neighbours of a phylogenetic tree in the rSPR graph by
#' the method described in Song (2003).
#'
#'
#' @param phy Phylogenetic tree.
#'
#' @return The number of neighbours of a phylogenetic tree in the rSPR Graph.
#' @export neighbour_size_song
#'
#' @examples
#' t <- ape::rtree(5, rooted = TRUE)
#' neighbour_size_song(t)
#' # Output:  [1] 26
#'
#' @references Song, Y. S. (2003). On the combinatorics of rooted binary
#' phylogenetic trees. \emph{ Annals of Com-binatorics, 7(3):365?37}
#'
#'
neighbour_size_song <- function(phy) {
  S_g <- compute_gamma(phy)
  n <- length(phy$tip.label)
  nei_size <- 2 * (n - 2) * (2 * n - 5) - 2 * S_g
  
  return(nei_size)
}

get_depths_and_sizes <- function(phy){
  n <- length(phy$tip.label)
  m <- phy$Nnode
  N <- dim(phy$edge)[1]
  phyr <- ape::reorder.phylo(phy, order = "postorder")
  phyr$edge.length <- rep(1, N)
  dmat <- ape::dist.nodes(x = phyr)
  descs <- sapply(1:(n + m),
                    function(i) {
                      phangorn::Descendants(x = phyr,
                                            node = i, type = "all")
                    }
    )
  raw.sizes <- unlist(lapply(descs, length))
  sizes <- raw.sizes + c(rep(0, n), rep(1, m))
  return(list(
    depths = dmat[n + 1, ],
    sizes = sizes
  ))
}

N_of_u <- function(phy, du, x){
  n <- ape::Ntip(phy)
  if(du<=0){
    ans <- 0
  }else{
    if(du == 1){
      ans <- 2*n - x - 3
    }else{
      ans <- 2*n - x - 5
    }
  }
}
#' Compute the neighbourhood size in the rSPR graph using the method from Matsen & Whidden (2017)
#'
#' @param phy a phylo object
#'
#' @return an integer counting how many neighbours a tree has in the rSPR graph.
#' @export neighbourhood_size_matsen
#'
neighbourhood_size_matsen <- function(phy){
  quantities <- get_depths_and_sizes(phy = phy)
  ds <- quantities$depths
  sizes <- quantities$sizes
  K <- length(sizes)
  Ns <- vector(length = K, mode = "numeric")
  for(k in 1:K){
    Ns[k] <- N_of_u(phy = phy,
                    du = ds[k], x = sizes[k])
  }
  return(sum(Ns))
}

#' Compute the size of rSPR neighbourhood
#'
#' @param phy Phylogenetic Tree
#' @param method Method used, could be "song" or "matsen"
#'  where \code{song} uses the method by Song (2003) and \code{matsen} uses the
#'  approach in Whidden and Matsen IV (2017). Default is method = "song".
#'
#' @return The size of the neighbourhood of a phylogenetic tree in the rSPR Graph
#' @export neighbourhood_size
#'
#' @examples
#' t <- ape::rtree(5, rooted = TRUE)
#' neighbourhood_size(t, method = "song")
#' neighbourhood_size(t, method = "matsen")
#' # Output:  [1] 150
#'
#'
#' @references Song, Y. S. (2003). On the combinatorics of rooted binary
#' phylogenetic trees. \emph{ Annals of Com-binatorics, 7(3):365?37}
#'
#' Whidden, C. and Matsen IV, F. A. (2017). Ricci-Ollivier
#' curvature of the rooted phylogenetic subtree-prune-regraft graph.
#' \emph{ Theoretical Computer Science, 699:1?2}
#'
neighbourhood_size <- function(phy, method = "song") {
  nei_size <- switch (method,
                      song = neighbour_size_song(phy),
                      matsen = neighbourhood_size_matsen(phy))
  
  return(nei_size)
}

#' Minimum rSPR neighbourhood size (i.e. for a ladder tree)
#'
#' @param n an integer giving the number of taxa.
#'
#' @return an integer counting the number of neighbours in the rSPR graph.
#' @export min_size
#'
min_size <- function(n){
  ans <- 3*n*n -13*n + 14
  return(ans)
}
min_size <- Vectorize(min_size)
#' Maximum rSPR neighbourhood size (i.e. for a balanced tree)
#'
#' @param n an integer giving the number of taxa.
#'
#' @return an integer counting the number of neighbours in the rSPR graph.
#' @export max_size

max_size <- function(n){
  S <- sum ( sapply(1:(n-2), function(m) floor(log2(m + 1))) )
  ans <- 4*(n-2)^2 - 2*S
  return(ans)
}
max_size <- Vectorize(max_size)

#' Numerically stable computation of the size of the neighbourhood of a ladder tree
#'
#' @param n an integer giving the number of taxa.
#' @param log logical. Whether to return the log of the requested quantity.
#'  Default is \code{FALSE}.
#'
#' @return an integer counting the size of the rSPR graph neighbourhood.
#' @export ladder_neighbourhood
ladder_neighbourhood <- function(n, log = FALSE){
  require(matrixStats)
  ## Lemma 5.1 in Whidden & Matsen (2017)
  if(n < 3) stop("n has to be greated than 3")
  require(matrixStats)
  l1 <- log(3) + 2*log(n)
  l2 <- log(13)  + log(n)
  ans <- log_diff_exp(matrixStats::logSumExp(c(l1, log(14))),
                      l2)
  if(!log) ans <- exp(ans)
  return(ans)
}
ladder_neighbourhood <- Vectorize(ladder_neighbourhood)

#' Numerically stable computation of the size of the neighbourhood of a balanced tree
#'
#' @param n an integer giving the number of taxa.
#' @param log logical. Whether to return the log of the requested quantity.
#'  Default is \code{FALSE}.
#'
#' @return an integer counting the size of the rSPR graph neighbourhood.
#' @export balanced_neighbourhood
balanced_neighbourhood <- function(n, log = FALSE){
  ## Lemma 5.1 in Whidden & Matsen (2017)
  if(n < 3) stop("n has to be >= 3")
  require(matrixStats)
  l1 <- log(4) + 2*log(n-2)
  lts <- sapply(1:(n-2), function(x) log(floor(log2(x + 1))))
  l2 <- log(2) + matrixStats::logSumExp(lts)
  ans <- log_diff_exp(l1, l2)
  if(!log) ans <- exp(ans)
  return(ans)
}
balanced_neighbourhood <- Vectorize(balanced_neighbourhood)
#' Numerically stable computation of the ratio of the size of the neighbourhood
#'  of a ladder tree to a balanced tree.
#'
#' @param n an integer giving the number of taxa.
#' @param log logical. Whether to return the log of the requested quantity.
#'  Default is \code{FALSE}.
#'
#' @return an integer counting the size of the rSPR graph neighbourhood.
#' @export Nsize_ratio
Nsize_ratio <- function(n, log = FALSE){
  ## Lemma 5.1 (i) of Whidden & Matsen (2017)
  if(n < 0) stop("n is supposed to be positive.")
  if (n <= 2) return(1)
  lnum <- log(3*n^2 - 13*n + 14)
  ldenom <- log(4*n^2 - 18*n + 20)   
  ans <- lnum-ldenom
  if(!log) ans <- exp(ans)
  return(ans)  
}
Nsize_ratio <- Vectorize(Nsize_ratio)