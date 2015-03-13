#' Compute a bootstrap standard error of the Harrell-Davis estimate of the qth quantile
#' 
#' @param x a numeric \code{vector}
#' @param q a desired quantile
#' @param nboot a number of bootstraps
#' @param SEED a logical indicating whether the seed should be set
#' @param cores a number of cores used for computation
#' @return the standard error of the Harrell-Davis estimate of the qth quantile
#' @export
hdseb <- function(x, q = .5, nboot = 100, SEED = TRUE, cores = 1) {
    if(SEED) {set.seed(2)}
    .hdseb(x, q, nboot, cores)
}