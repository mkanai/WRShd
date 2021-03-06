#' Compute the Harrell-Davis estimate of the qth quantile
#' 
#' @param x a numeric \code{vector}
#' @param q a desired quantile
#' @param na.rm a logical indicating whether missing values should be checked for removal
#' @param cores a number of cores used for computation
#' @return the Harrell-Davis estimate of the qth quantile of \code{x}
#' @export
hd <- function(x, q = .5, na.rm = TRUE, cores = 1) {
    .hd(x, q, na.rm, cores)
}