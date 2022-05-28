#' Excise a clade from a tree
#'
#' @param phy a phylo object.
#' @param clade a list of tip names in a clade.
#'
#' @return a list with two phylo objects: one is the tree without the clade and
#'  the other is the subtree with the clade.
#' @export f_c
#'
f_c <- function(phy, clade){
  ## x prime
  wo.clade <- ape::drop.tip(phy = phy, tip = clade,
                       trim.internal = FALSE,
                       subtree = TRUE, rooted = TRUE)
  pos <- grep("_tip", wo.clade$tip.label)
  if(length(pos) > 1) stop("Taxa in 'clade' are not monophyletic")
  wo.clade$tip.label[pos] <- "l"
  ## remaining subtree
  other.side <- ape::drop.tip(phy = phy,
                         tip = setdiff(phy$tip.label, clade),
                         trim.internal = TRUE,
                         subtree = FALSE, rooted = TRUE)
  return(
    list(
      x_prime = wo.clade,
      phi_x_c = other.side
    )
  )
}

#' Safe version of \code{f_c}
#'
#' @param phy a phylo object.
#' @param clade a list of tip names in a clade.
#'
#' @return a list with two phylo objects: one is the tree without the clade and
#'  the other is the subtree with the clade. Or \code{NA} if
#'   the clade is not compatible with \code{phy}
#' @export safe_f_c
#'
safe_f_c <- function(phy, clade){
  tryCatch(
    f_c(phy, clade),
    error = function(e){
      return(NA)
    }
  )
}

#' Get the largest clade in a tree
#'
#' @param phy a phylo object.
#'
#' @return a list containing the largest tree and its size
#' @export get_largest_clade
#'
get_largest_clade <- function(phy){
  
  nt <- ape::Ntip(phy)
  nodes <- 1:ape::Nedge(phy)
  nodes <- nodes[nodes > nt + 1]
  clades <-lapply(nodes, function(i){
    phangorn::Descendants(x = phy, node = i, type = "tips")  
  })
  
  clade.sizes <- unlist(lapply(clades, function(x) length(x[[1]]))) 
  
  pos <- nodes[which.max(clade.sizes)]
  
  largest.clade <- ape::extract.clade(phy = phy,
                                 node = pos)$tip.label
  list(
    clade = largest.clade,
    size = length(largest.clade)
  )
}


#' Get a random clade
#'
#' @param phy a phylo object
#'
#' @return a list containing a random clade and its size.
#' @export get_random_clade
get_random_clade <- function(phy){
  nt <- ape::Ntip(phy)
  nodes <- 1:ape::Nedge(phy)
  nodes <- nodes[nodes > nt + 1]
  pos <-  sample(nodes, 1, replace = FALSE)
  r.clade <- ape::extract.clade(phy = phy,
                           node = pos)$tip.label
  return(
    list(
      clade = r.clade,
      size = length(r.clade)
    )
  )
}