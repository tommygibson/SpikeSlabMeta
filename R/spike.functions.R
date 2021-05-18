### Functions to be used in a bunch of things


# 3 significant digits (won't round to 2)
sigfig <- function(x, n=3){ 
  
  trimws(format(round(x, n), nsmall = n), which = "both")
  
}   

# summary from MCMC sims for CTSs
make_CTS_sum1 <- function(x, M = 100){
  
  # averaging over stat_new in each iteration to get the posterior
  stats <- apply(x[, 1:M], 1, mean)
  
  
  means <- mean(stats)
  low.up <- quantile(stats, c(.025, .5, .975))
  sds <- sd(stats)
  stat.summ <- c(means, sds, low.up)
  
  return(stat.summ)
  
}

CTS.overall.sum <- function(x, target){
  
  # LR+
  bias <- mean(x[,1]) - target[1] # bias
  var.est <- var(x[,1]) # mean of estimator
  avg.sd <- mean(x[,2]) # average SD
  cover.95 <- sum(x[,3] < target[1] & x[,5] > target[1]) / dim(x)[1] # 95% coverage probability
  length.95 <- mean(x[, 5] - x[, 3]) # average 95% length
  RMSE <- sqrt(mean((x[, 1] - target[1])^2)) # root(MSE)
  
  return(c(bias, var.est, avg.sd, cover.95, length.95, RMSE))
}

reject_by_cutoff <- function(cut, mat){
  
  return(apply(mat, 2, function(x) sum(x > cut) / length(x)))
  
}

#### Init generators

init.gen <- function(){
  list(
    delta = runif(S, -1, 1),
    beta = runif(S, -1, 1),
    nu = runif(S, -1, 1),
    delta0 = runif(1, -1, 1),
    beta0 = runif(1, -1, 1),
    nu0 = runif(1, -1, 1),
    sigma.delta = runif(1, 0.2, 1),
    sigma.beta = runif(1, 0.2, 1),
    sigma.nu = runif(1, 0.2, 1),
    deltanew = runif(M, -1, 1),
    betanew = runif(M, -1, 1),
    nunew = runif(M, -1, 1)
  )
}

init.gen.spike <- function(){
  list(
    delta = runif(S, -1, 1),
    beta = runif(S, -1, 1),
    nu = runif(S, -1, 1),
    delta1 = runif(1, -1, 1),
    beta0 = runif(1, -1, 1),
    nu0 = runif(1, -1, 1),
    sigma.delta = runif(1, 0.2, 1),
    sigma.beta = runif(1, 0.2, 1),
    sigma.nu = runif(1, 0.2, 1),
    deltanew = runif(M, -1, 1),
    betanew = runif(M, -1, 1),
    nunew = runif(M, -1, 1),
    rho = 1
  )
}