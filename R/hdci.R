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
    x <- .na.omit(x)
    if(pr && sum(duplicated(x) > 0)) {print("Duplicate values detected; use hdpb")}
    se <- hdseb(x, q, nboot, cores)
    crit <- .5064 / (length(x)^(.25)) + 1.96
    if ((q<=.2 || q>=.8) && length(x) <= 20) {
        crit <- (-6.23) / length(x) + 5.01
    }
    if((q<=.1 || q>=.9) && length(x) <= 40) {
        crit <- 36.2 / length(x) + 1.31
    }
    if(length(x)<=10){
        print("The number of observations is less than 11.")
        print("Accurate critical values have not been determined for this case.")
    }
    low <- hd(x, q, cores) - crit*se
    hi <- hd(x, q, cores) + crit*se
    list(ci = c(low, hi), crit = crit, se = se)
}

