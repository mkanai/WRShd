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
#' @export

qcomhd <- function(x, y, q = c(.1, .25, .5, .75, .9), nboot = 2000, plotit = TRUE, SEED = TRUE, xlab = "Group 1", ylab = "Est.1-Est.2", alpha = .05, cores = 1) {
    if(SEED) {set.seed(2)}
    pv <- NULL
    output <- matrix(NA, nrow = length(q), ncol = 10)
    dimnames(output) <- list(NULL, c("q", "n1", "n2", "est.1", "est.2", "est.1_minus_est.2", "ci.low", "ci.up", "p_crit", "p-value"))
    for (i in 1:length(q)) {
        output[i,1] <- q[i]
        output[i,2] <- length(.na.omit(x))
        output[i,3] <- length(.na.omit(y))
        output[i,4] <- hd(x, q = q[i])
        output[i,5] <- hd(y, q = q[i])
        output[i,6] <- output[i,4] - output[i,5]
        temp <- .qcom.sub(x, y, nboot = nboot, q = q[i], alpha = alpha, cores = cores)
        output[i,7] <- temp$ci[1]
        output[i,8] <- temp$ci[2]
        output[i,10] <- temp$p.value                                                                                                      
    }                                                                                                                              
    temp <- order(output[,10], decreasing = TRUE)                                                                                        
    zvec <- alpha / c(1:length(q))                                                                                                      
    output[temp,9] <- zvec                                                                                                            
    output <- data.frame(output)                                                                                                   
    output$signif <- rep("YES", nrow(output))
    for(i in 1:nrow(output)) {
        if(output[temp[i],10] > output[temp[i],9])output$signif[temp[i]]="NO"                                                            
        if(output[temp[i],10] <= output[temp[i],9])break                                                                                 
    }                                                                                                                              
    if(plotit) {                                                                                                                    
        xax <- rep(output[,4],3)                                                                                                          
        yax <- c(output[,6], output[,7], output[,8])                                                                                        
        plot(xax, yax, xlab = xlab, ylab = ylab, type="n")                                                                                     
        points(output[,4], output[,6], pch = "*")                                                                                          
        lines(output[,4], output[,6])                                                                                                   
        points(output[,4], output[,7], pch = "+")                                                                                          
        points(output[,4], output[,8], pch = "+")                                                                                          
    }                                                                                                                              
    output                                                                                        
}                     
