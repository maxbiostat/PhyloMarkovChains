#' Create a symmetric matrix from a vector
#'
#' @param x a vector
#' @param diag boolean indicating whether to include the diagonal (default is \code{TRUE})
#' @param byrow boolean indicating whether to compute the matrix by row.
#'
#' @return a symmetric matrix
#' @export vec2symMat
#'
vec2symMat <- function(x, diag = TRUE, byrow = FALSE) {
  # stolen from https://github.com/mikewlcheung/metasem/blob/17542243fc3f69c1b6d3c53ec68c00f4f8bbb81f/R/vec2symMat.R
  m <- length(x)
  d <- if (diag)  1 else - 1
  n <- floor((sqrt(1 + 8 * m) - d) / 2)
  if (m != n * (n + d) / 2)
    stop("Cannot make a square matrix as the length of \"x\" is incorrect.")
  mat <- diag(n)
  
  ## Row major
  if (byrow) {
    mat[upper.tri(mat, diag = diag)] <- x
    index <- lower.tri(mat)
    mat[index] <- t(mat)[index]
  } else {
    ## Column major: default behavior
    mat[lower.tri(mat, diag = diag)] <- x
    # Just mirroring the matrix, exclude the diagonals
    ## mat[upper.tri(mat, diag=FALSE)] <- mat[lower.tri(mat, diag=FALSE)]
    ## Corrected a bug
    index <- upper.tri(mat)
    mat[index] <- t(mat)[index]
  }
  mat
}

#' Get the shapes of all trees in a set
#'
#' @param all_trees the set of trees to be compared.
#'  Usually the output of phangorn::allTrees()
#' @param ncores number of cores to be used. Default is \code{ncores=6}.
#'
#' @return a binary matrix which is 1 if two trees have the same shape and 0 otherwise
#' @export get_shapes
#'
get_shapes <- function(all_trees, ncores = 6) {
  
  requireNamespace(foreach)
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  K <- length(all_trees)
  grid <- subset(expand.grid(1:K, 1:K), Var1 < Var2)
  N <- nrow(grid)
  comps <- foreach::foreach(k = 1:N, 
                            .combine = "cbind") %dopar% {
                              i <- grid[k, 1]
                              j <- grid[k, 2]
                              ape::all.equal.phylo(all_trees[[i]],
                                                   all_trees[[j]],
                                                   use.edge.length = FALSE,
                                                   use.tip.label = FALSE)
                            }
  
  
  parallel::stopCluster(cl)
  
  comps <- as.numeric(comps)
  Shape_num <- vec2symMat(x = comps, diag = FALSE, byrow = TRUE)
  diag(Shape_num) <- 1
  return(Shape_num)
}

##############

#' Project (lump) a transition matrix on the space of trees to the space of shapes
#'
#' @param transMat Markov transtion matrix on the space of trees
#' @param n the number of taxa.
#' @param ncores the number of cores to be used.
#'
#' @return a projected (lumped) transition matrix on the space of shapes
#' @export shape_lump
#'
shape_lump <- function(transMat, n, ncores = 6) {
  
  all.trees <- phangorn::allTrees(n, rooted = TRUE)
  K <- length(all.trees)
  Shape_num <- get_shapes(all.trees, ncores = ncores)
  
  dif_shapes <- unique.array(Shape_num)
  DS <- nrow(dif_shapes)
  shapes_each_total <- rowSums(dif_shapes)
  Max_shape <- max(shapes_each_total)
  Mat_shapes <- matrix(0, nrow = Max_shape, ncol = DS)
  
  
  for (i in 1:DS) {
    Mat_shapes[(1:shapes_each_total[i]), i] <- which(apply(Shape_num, 1,
                                                           function(x)
                                                             all.equal(x[],
                                                                       dif_shapes[i,]) == TRUE))
    
  }
  
  one_shape_each <- array(NA, dim = DS)
  
  for (i in 1:DS) {
    one_shape_each[i] <- Mat_shapes[1, i]
  }
  
  Mat_lump <- matrix(NA, nrow = DS, ncol = DS)
  
  for (i in 1:DS) {
    x_S <- array(NA, dim = K)
    x_S <- transMat[one_shape_each[i], ]
    for (j in 1:DS) {
      Mat_lump[i, j] <- sum(x_S[Mat_shapes[, j]])
    }
  }
  
  return(Mat_lump)
}