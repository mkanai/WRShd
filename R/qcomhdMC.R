#' Compare quantiles via the Harrell-Davis estimator
#' 
#' @param x a numeric \code{vector}
#' @param y a numeric \code{vector}
#' @param q a desired quantile
#' @param nboot a number of bootstraps
#' @param plotit 
#' @param SEED
#' @param xlab
#' @param ylab
#' @param alpha a significance level
#' @param cores a number of cores used for computation
#' @export
qcomhdMC <- function(x, y, q = c(.1, .25, .5, .75, .9), nboot = 2000, plotit = TRUE, SEED = TRUE, xlab = "Group 1", ylab = "Est.1-Est.2", alpha = .05, cores = -1) {
    require(parallel)
    if(cores <= 0){cores <- parallel::detectCores()}
    qcomhd(x, y, q, nboot, plotit, SEED, xlab, ylab, alpha, cores)
}