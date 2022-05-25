#' Log of 1-x
#'
#' @param x a vector of values between 0 and 1 for which to compute log(1 - x).
#'
#' @return log(1 - x).
#' @export log1m
#'
log1m <- function(x) {
  return(log1p(-x))
}
#' Log of 1-exp(x)
#'
#' @param a a vector of values a such that exp(a) < 1
#'
#' @return log(1-exp(a))
#' @export 
#'
log1m_exp <- function(a) {
  if (a > -0.693147) {
    return(log(-expm1(a))) 
  } else {
    return(log1m(exp(a)));
  }
}
#' Log of exp(x)-exp(y)
#'
#' @param x a scalar
#' @param y a scalar
#'
#' @return log(exp(x)-exp(y))
#' @export log_diff_exp
#'
log_diff_exp <- function (x, y){
  return(x + log1m_exp(y - x))
}