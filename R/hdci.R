#' Compute a 1-\code{alpha} confidence for the Harrell-Davis estimate of the qth quantile
#' 
#' @param x a numeric \code{vector}
#' @param q a desired quantile
#' @param alpha a significance level
#' @param nboot a number of bootstraps
#' @param SEED a logical indicating whether the seed should be set
#' @param pr a logical indicating whether
#' @param cores a number of cores used for computation
#' @return the 1-\code{alpha} confidence for the Harrell-Davis estimate of the qth quantile
#' @export

hdci <- function(x, q = .5, alpha = .05, nboot = 100, SEED = TRUE, pr = TRUE, cores = 1) {
    if(SEED)set.seed(2)
    if(alpha != .05) {stop("Use the function qcipb. Generally works well even when alpha is not equal to .05")}
    .hdci(x, q, nboot, pr, cores)
}

