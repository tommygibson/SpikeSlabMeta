### Functions to be used in a bunch of things


# 3 significant digits (won't round to 2)
sigfig <- function(x, n = 3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
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
  var.est <- var(x[,1]) # variance of estimator
  avg.sd <- mean(x[,2]) # average SD
  cover.95 <- sum(x[,3] < target[1] & x[,5] > target[1]) / dim(x)[1] # 95% coverage probability
  length.95 <- mean(x[, 5] - x[, 3]) # average 95% length
  RMSE <- sqrt(mean((x[, 1] - target[1])^2)) # root(MSE)
  
  return(c(bias, var.est, avg.sd, cover.95, length.95, RMSE))
}

reject_by_cutoff <- function(cut, mat){
  
  return(apply(mat, 2, function(x) sum(x > cut) / length(x)))
  
}

### functions for calculating CTSs from hyperparameter estimates

cts0_from_gamma <- function(gamma, M){
  beta0 <- gamma[1]
  delta0 <- gamma[2]
  nu0 <- gamma[3]
  sigma.beta <- gamma[4]
  sigma.delta <- gamma[5]
  sigma.nu <- gamma[6]
  
  beta.new <- rnorm(M, mean = beta0, sd = sigma.beta)
  delta.new <- rnorm(M, mean = delta0, sd = sigma.delta)
  nu.new <- rnorm(M, mean = nu0, sd = sigma.nu)
  
  pi1.new <- 1 / (1 + exp(-(beta.new + delta.new / 2)))
  pi0.new <- 1 / (1 + exp(-(beta.new - delta.new / 2)))
  psi.new <- 1 / (1 + exp(-(nu.new)))
  
  sens.new <- (pi1.new * psi.new) / (pi1.new * psi.new + pi0.new * ( 1 - psi.new))
  spec.new <- ((1 - pi0.new) * (1 - psi.new)) / ((1 - pi0.new) * (1 - psi.new) + (1 - pi1.new) * psi.new)
  LRm.new <- (1 - sens.new) / spec.new
  LRp.new <- sens.new / (1 - spec.new)
  
  mc.est <- apply(cbind(LRm.new, LRp.new, 1 - pi0.new, pi1.new, sens.new, spec.new), MARGIN = 2, mean)
  return(mc.est)
}

ctsh_from_gamma <- function(gamma){
  beta0 <- gamma[1]
  delta0 <- gamma[2]
  nu0 <- gamma[3]
  sigma.beta <- gamma[4]
  sigma.delta <- gamma[5]
  sigma.nu <- gamma[6]
  
  pi1 <- 1 / (1 + exp(-(beta0 + delta0 / 2)))
  pi0 <- 1 / (1 + exp(-(beta0 - delta0 / 2)))
  psi <- 1 / (1 + exp(-nu0))
  
  sens <- pi1 * psi / (pi1 * psi + pi0 * (1 - psi))
  spec <- (1 - pi0) * (1 - psi) / ((1 - pi0) * (1 - psi) + (1 - pi1) * psi)
  
  LRp <- sens / (1 - spec)
  LRm <- (1 - sens) / spec
  
  hyper.est <- c(LRm, LRp, 1 - pi0, pi1, sens, spec)
  
  return(hyper.est)
  
}

make_cts0_from_hyper <- function(hyperparameter_matrix, M){
  t(apply(hyperparameter_matrix, MARGIN = 1, cts0_from_gamma, M = M))
}

make_ctsh_from_hyper <- function(hyperparameter_matrix){
  t(apply(hyperparameter_matrix, MARGIN = 1, ctsh_from_gamma))
}

simple_summary <- function(vec){
  return(c(mean(vec), sd(vec), quantile(vec, c(.025, .5, .975))))
}

cts_summary_from_gamma <- function(hyperparameter_matrix, M){

  summary_0 <- t(apply(make_cts0_from_hyper(hyperparameter_matrix, M = M), MARGIN = 2, simple_summary))
  summary_h <- t(apply(make_ctsh_from_hyper(hyperparameter_matrix), MARGIN = 2, simple_summary))
  
  rownames(summary_0) <- rownames(summary_h) <- NULL
  
  return(cbind(rep(c("LRm", "LRp", "NPV", "PPV", "sens", "spec"), 2),
               rep(c("CTS0", "CTSh"), each = 6),
               rbind(summary_0, summary_h)))
  
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