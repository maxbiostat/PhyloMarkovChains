#' Run a (phylogenetic) Markov chain from a transition matrix
#'
#' @param Niter Integer number of iterations
#' @param MH_mat Transition matrix
#' @param ini.state initial state. If \code{NULL}, a state is chosen uniformly.
#'
#' @return A markov chain of tree indices
#' @export run_one_MC
#'
run_one_MC <- function(Niter, MH_mat, ini.state = NULL){
  
  K <- nrow(MH_mat)
  if(ncol(MH_mat) != K) stop("MH_mat is not square")
  
  States <- paste0("tree_", 1:K)
  
  phylo.MH <- methods::new( ## needed for markovchainSequence()
    "markovchain",
    states = States,
    transitionMatrix = MH_mat,
    name = "rootedMH"
  )
  
  if(is.null(ini.state)) ini.state <- sample(States, 1)
  
  out <- markovchain::markovchainSequence(n = Niter,
                                          markovchain = phylo.MH,
                                          t0 = ini.state)
  return(out)
}