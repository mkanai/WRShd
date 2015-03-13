#' Compute a bootstrap 1-\code{alpha} confidence for the Harrell-Davis estimate of the qth quantile
#' 
#' @param x a numeric \code{vector}
#' @param q a desired quantile
#' @param alpha a significance level
#' @param nboot a number of bootstraps
#' @param SEED a logical indicating whether the seed should be set
#' @param nv a null-value when computing a p-value
#' @param cores a number of cores used for computation
#' @return the a bootstrap 1-\code{alpha} confidence for the Harrell-Davis estimate of the qth quantile
#' @export
hdpb <- function(x, q = .5, alpha = .05, nboot = 2000, SEED = TRUE, nv = 0, cores = 1) {
    if(SEED) {set.seed(2)}
    .hdpb(x, q, alpha, nboot, nv, cores)
}
