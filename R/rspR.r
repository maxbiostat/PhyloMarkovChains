############
### Wrappers for the great rpsr (https://github.com/cwhidden/rspr) by Chris Whidden (FHCRC)
### Note that the functions here assume rspr is installed in the system, ie., at /usr/bin/rspr
############
#' Compute the rooted SPR (rSPR) distance between two trees
#'
#' @param t1 a phylo object.
#' @param t2 a phylo object.
#' @param exact logical. Whether to compute the exact distance (default is \code{TRUE}) or 
#' an approximation.
#'
#' @return a numerical distance
#' @export rspr
#'
#' @details  Note that this assumes rspr, \url{https://github.com/cwhidden/rspr}
#'  is installed in the system, ie., at /usr/bin/rspr.
rspr <- function(t1, t2, exact = TRUE){
  tts <- list(t1, t2)
  class(tts) <- "multiPhylo"
  raw <- system2(input = paste(ape::write.tree(tts)), command = "rspr", stdout = TRUE)
  if(exact){    
    d <- strsplit(raw[grep("total exact", raw)], "=")[[1]][2]
  }else{
    d <- strsplit(raw[grep("approx drSPR", raw)], "=")[[1]][2]
  }
  return(as.numeric(d))
}

#' Compute a matrix of rSPR distances
#'
#' @param trees a multiPhylo object
#' @param type if 'restricted', this will compute distances up to \code{maxdist}.
#' If \code{onerow}, will compute the set of SPR distances from the first tree and
#' if \code{full} ompute the full matrix (type = "full").
#' While fun for small data sets, may be impractical for large ones. 
#' @param maxdist numerical. The maximum distance up to which to compute distances.
#'
#' @return a matrix of distances
#' @export rspr_matrix
#'
rspr_matrix <- function(trees,
                        type = c("restricted", "onerow", "full"),
                        maxdist = 1){ 
  type <- match.arg(type)
  
  if(!inherits(trees, "multiPhylo")) stop("Input is not of class multiPhylo") 
  
  root.test <- unlist(lapply(trees, ape::is.rooted)) # testing if trees are rooted
  tmp <- tempfile("rsprmatrix", fileext = ".csv") # Unfortunately, we'll have to resort to this hack for the time being.
  if(type == "restricted"){
    if(any(root.test)){ # If any of the trees in the set is rooted, uses rooted (default) version
      dists <- system2(input = paste(ape::write.tree(trees)), command = "rspr",
                       args = paste("-pairwise_max", maxdist), stdout = tmp)
    }else{
      dists <- system2(input = paste(ape::write.tree(trees)), command = "rspr",
                       args = paste("-unrooted -pairwise_max", maxdist),  stdout = tmp)
    }    
  }else{
    if(type == "onerow"){
      if(any(root.test)){ 
        dists <- system2(input = paste(ape::write.tree(trees)),
                         command = "rspr",  args = "-pairwise 0 1", stdout = tmp)
      }else{
        dists <- system2(input = paste(ape::write.tree(trees)),
                         command = "rspr",  args = "-unrooted -pairwise 0 1", stdout = tmp)
      }
    }else{
      if(any(root.test)){
        dists <- system2(input = paste(ape::write.tree(trees)),
                         command = "rspr",  args = "-pairwise", stdout = tmp)
      }else{
        dists <- system2(input = paste(ape::write.tree(trees)),
                         command = "rspr",  args = "-unrooted -pairwise", stdout = tmp)
      }
    }
  }
  m <- as.matrix(utils::read.table(tmp, sep = ",", header = FALSE, fill = TRUE))
  if(any(dim(m)>1))  m[lower.tri(m)] <- t(m)[lower.tri(t(m))]
  rownames(m) <- colnames(m) <- NULL
  return(m)
}