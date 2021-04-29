######### simulation for spike/slab meta-analysis

set.seed(1234)

library(R2jags)
library(knitr)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(here)

# only necessary if hasn't been imported before
# font_import(pattern = "lmroman*")

#### First, the models

####### STANDARD, NO SPIKE/SLAB
sink("meta_confusion.txt")
cat("
model
{
  for(i in 1:S){
    y[i,1] ~ dbin(pi[i,1], n[i,1])
    y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    # meta-regression on the log-odds ratio
    # do we include intercept?
    
    delta[i] ~ dnorm(delta0, D.delta)
    
    nu[i] ~ dnorm(nu0, D.nu)
    
    beta[i] ~ dnorm(beta0, D.beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  delta0 ~ dnorm(c, d)
  nu0 ~ dnorm(e, f)
  
  # HALF T ON ALL THESE BITCHES
  D.beta <- pow(sigma.beta, -2)
  D.nu <- pow(sigma.nu, -2)
  D.delta <- pow(sigma.delta, -2)
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 1) T(0,)
  sigma.nu ~ dt(0, 1, 1) T(0,)
  sigma.delta ~ dt(0, 1, 1) T(0,)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D.delta)
    betanew[j] ~ dnorm(beta0, D.beta)
    nunew[j] ~ dnorm(nu0, D.nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi1new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))
    pi0new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))
    
    pi11new[j] <- pi1new[j] * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - pi1new[j]) * psinew[j]            # P(no event, risk)
    pi01new[j] <- pi0new[j] * (1 - psinew[j])            # P(event, no risk)
    pi00new[j] <- (1 - pi0new[j]) * (1 - psinew[j])      # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
  pi1.h <- 1 / (1 + exp(-(beta0 + delta0 / 2)))
  pi0.h <- 1 / (1 + exp(-(beta0 - delta0 / 2)))
  psi.h <- 1 / (1 + exp(-(nu0)))
  
  pi11.h <- pi1.h * psi.h
  pi10.h <- (1 - pi1.h) * psi.h
  pi01.h <- pi0.h * (1 - psi.h)
  pi00.h <- (1 - pi0.h) * (1 - psi.h)
  
  sens.h <- pi11.h / (pi11.h + pi01.h)
  spec.h <- pi00.h / (pi00.h + pi10.h)
  
  LRp.h <- sens.h / (1 - spec.h)
  LRm.h <- (1 - sens.h) / spec.h
  
}", fill = TRUE)
sink()

########## SPIKE/SLAB
sink("meta_confusion_spike.txt")
cat("
model
{
  for(i in 1:S){
    y[i,1] ~ dbin(pi[i,1], n[i,1])
    y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    
    delta[i] ~ dnorm(delta0, D.delta)
    
    nu[i] ~ dnorm(nu0, D.nu)
    
    beta[i] ~ dnorm(beta0, D.beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  nu0 ~ dnorm(e, f)
  
  delta0 <- delta1 * rho
  delta1 ~ dnorm(c, d)
  rho ~ dbern(p)
  spike <- 1 - rho
  
  # HALF T ON ALL THESE BITCHES
  D.beta <- pow(sigma.beta, -2)
  D.nu <- pow(sigma.nu, -2)
  D.delta <- pow(sigma.delta, -2)
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 1) T(0,)
  sigma.nu ~ dt(0, 1, 1) T(0,)
  sigma.delta ~ dt(0, 1, 1) T(0,)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D.delta)
    betanew[j] ~ dnorm(beta0, D.beta)
    nunew[j] ~ dnorm(nu0, D.nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi1new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))
    pi0new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))
    
    pi11new[j] <- pi1new[j] * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - pi1new[j]) * psinew[j]            # P(no event, risk)
    pi01new[j] <- pi0new[j] * (1 - psinew[j])            # P(event, no risk)
    pi00new[j] <- (1 - pi0new[j]) * (1 - psinew[j])      # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
}", fill = TRUE)
sink()

## some useful functions

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
########### DATA GENERATION

# FIND EXPECTED VALUES FOR STATS GIVEN HYPERPARAMETER VALUES

nsim <- 10000
sigma.nu <- sigma.beta <- 0.4
nu0 <- beta0 <- log(0.15 / (1 - 0.15))
delta0 <- 2
delta0fix <- 0
sigma.delta <- 0.4

deltasims <- rnorm(nsim, mean = delta0, sd = sigma.delta)
nusims <- rnorm(nsim, mean = nu0, sd = sigma.nu)
betasims <- rnorm(nsim, mean = beta0, sd = sigma.beta)


pi11sims <- 1 / (1 + exp(-(betasims + deltasims / 2))) * 1 / (1 + exp(-nusims))
pi10sims <- (1 - 1 / (1 + exp(-(betasims + deltasims / 2)))) * 1 / (1 + exp(-nusims))
pi01sims <- 1 / (1 + exp(-(betasims - deltasims / 2))) * (1 - 1 / (1 + exp(-nusims)))
pi00sims <- (1 - 1 / (1 + exp(-(betasims - deltasims / 2)))) * (1 - 1 / (1 + exp(-nusims)))

pi11fix <- 1 / (1 + exp(-(betasims + delta0fix / 2))) * 1 / (1 + exp(-nusims))
pi10fix <- (1 - 1 / (1 + exp(-(betasims + delta0fix / 2)))) * 1 / (1 + exp(-nusims))
pi01fix <- 1 / (1 + exp(-(betasims - delta0fix / 2))) * (1 - 1 / (1 + exp(-nusims)))
pi00fix <- (1 - 1 / (1 + exp(-(betasims - delta0fix / 2)))) * (1 - 1 / (1 + exp(-nusims)))

senssims <- pi11sims / (pi11sims + pi01sims)
specsims <- pi00sims / (pi00sims + pi10sims)

LRpsims <- senssims / (1 - specsims)
LRmsims <- (1 - senssims) / specsims

PPVsims <- pi11sims / (pi11sims + pi10sims)
NPVsims <- pi00sims / (pi00sims + pi01sims)

sensfix <- pi11fix / (pi11fix + pi01fix)
specfix <- pi00fix / (pi00fix + pi10fix)

LRpfix <- sensfix / (1 - specfix)
LRmfix <- (1 - sensfix) / specfix

PPVfix <- pi11fix / (pi11fix + pi10fix)
NPVfix <- pi00fix / (pi00fix + pi01fix)

RDfix <- PPVfix - (1 - NPVfix)
RRfix <- PPVfix / (1 - NPVfix)

stats.sims <- cbind.data.frame(LRmsims, LRpsims, NPVsims, PPVsims, senssims, specsims)
stats.fix <- cbind.data.frame(LRmfix, LRpfix, NPVfix, PPVfix, sensfix, specfix)
names(stats.sims) <- c("LRm", "LRp", "NPV", "PPV", "Sens", "Spec")
names(stats.fix) <- c("LRm", "LRp", "NPV", "PPV", "Sens", "Spec")

apply(stats.sims, 2, mean)
apply(stats.sims, 2, quantile, c(.025, .5, .975))

apply(stats.fix, 2, mean)
apply(stats.fix, 2, quantile, c(.025, .5, .975))

stats.means <- apply(stats.sims, 2, mean)


###### Okay, now actual simulation
###### need to generate 2x2 tables and then apply models to the data repeatedly


###### Spike/slab model
###### Assess how well it catches true zeros and avoids real nonzeros

#### FIRST, HOW WELL IT CATCHES ZEROS

S <- 10
N <- 500
K <- 1000
sigma.delta <- c(0.1, 0.25, 0.5)

spike.summary <- matrix(nrow = K, ncol = length(sigma.delta))
set.seed(10101)
for(i in 1:length(sigma.delta)){
  for(k in 1:K){
    
    delta.i <- rnorm(S, mean = 0, sd = sigma.delta[i])
    beta.i <- rnorm(S, mean = beta0, sd = sigma.beta)
    nu.i <- rnorm(S, mean = nu0, sd = sigma.nu)
    
    pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
    pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
    pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
    pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))
    
    probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities
    
    tables <- t(apply(probs, 1, rmultinom, n = 1, size = N)) # generate contingency tables
    
    # prepare data for going into jags
    y <- tables[,c(1, 3)]
    n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
    n.tot <- rep(N, S)
    
    M <- 100
    
    meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                     a = -1.5, b = 0.5, c = 0, d = 0.5, e = -1.5, f = 0.5, p = 0.5)
    
    
    # as a first run we'll just follow the hyperparameters
    meta.params <- c("spike")
    
    spike.zero <- jags(data = meta.dat, inits = init.gen.spike, parameters.to.save = meta.params,
                       model.file = "meta_confusion_spike.txt",
                       n.chains = 2, n.iter = 6000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
    
    # columns 1 and 3 are mean (1 - rho)
    # columns 2 and 4 are indicators for mean (1 - rho) < x
    spike.summary[k, i] <- spike.zero$BUGSoutput$summary[1]
    
    
  }
}

seed.after.zeros <- .Random.seed

####### SECOND: HOW WELL IT CATCHES NON-ZEROS


K <- 1000
sigma.delta <- c(0.1, 0.25, 0.5)
delta0.nonzero <- c(1, 2)

spike.summary.nonzero <- matrix(nrow = K, ncol = length(sigma.delta) * length(delta0.nonzero))

.Random.seed <- seed.after.zeros

for(i in 1:length(sigma.delta)){
  for(j in 1:length(delta0.nonzero)){
    for(k in 1:K){
      
      delta.i <- rnorm(S, mean = delta0.nonzero[j], sd = sigma.delta[i])
      beta.i <- rnorm(S, mean = beta0, sd = sigma.beta)
      nu.i <- rnorm(S, mean = nu0, sd = sigma.nu)
      
      pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
      pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
      pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
      pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))
      
      probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities
      
      tables <- t(apply(probs, 1, rmultinom, n = 1, size = N)) # generate contingency tables
      
      # prepare data for going into jags
      y <- tables[,c(1, 3)]
      n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
      n.tot <- rep(N, S)
      M <- 100
      
      meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                       a = -2, b = 0.5, c = 0, d = 0.5, e = -1, f = 0.5, p = 0.5)
      
      
      # as a first run we'll just follow the hyperparameters
      meta.params <- c("spike")
      
      spike.zero <- jags(data = meta.dat, inits = init.gen.spike, parameters.to.save = meta.params,
                         model.file = "meta_confusion_spike.txt",
                         n.chains = 2, n.iter = 6000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
      
      # columns represent the two values for sigma_delta (small and moderate)
      
      spike.summary.nonzero[k, (length(delta0.nonzero) * (i - 1) + j)] <- spike.zero$BUGSoutput$summary[1]
      
      
    }
  }
}

seed.after.nonzero <- .Random.seed



####### THIRD: BIAS, COVERAGE %, 95% LENGTH FOR CTSs

# K determined by MCSE of bias for LR+ 
# Var(\hat{\theta}) \le 0.25, want MCSE < 0.01
K <- 2500

CTSs <- list()
for(i in 1:6){
  CTSs[[i]] <- matrix(nrow = K, ncol = 5)
}

names(CTSs) <- c("LRm", "LRp", "NPV", "PPV", "sens", "spec")

delta0 <- 2
sigma.delta <- 0.4

.Random.seed <- seed.after.nonzero
for(k in 1:K){
    
  delta.i <- rnorm(S, mean = delta0, sd = sigma.delta)
  beta.i <- rnorm(S, mean = beta0, sd = sigma.beta)
  nu.i <- rnorm(S, mean = nu0, sd = sigma.nu)
  
  pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
  pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
  pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
  pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))
  
  probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities
  
  tables <- t(apply(probs, 1, rmultinom, n = 1, size = N)) # generate contingency tables
  
  # prepare data for going into jags
  y <- tables[,c(1, 3)]
  n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
  n.tot <- rep(N, S)
  M <- 100
  
  meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                   a = -1, b = 0.5, c = 0, d = 0.5, e = -1, f = 0.5)

  meta.params <- c("LRmnew", "LRpnew", "sensnew", "specnew", "PPVnew", "NPVnew")
  
  std.model <- jags(data = meta.dat, inits = init.gen, parameters.to.save = meta.params,
                    model.file = "meta_confusion.txt",
                    n.chains = 2, n.iter = 6000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  
  
  CTSs[[1]][k,] <- make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,1:M], M)
  CTSs[[2]][k,] <- make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(M + 1):(2 * M)], M)
  CTSs[[3]][k,] <- make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(2 * M + 1):(3 * M)], M)
  CTSs[[4]][k,] <- make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(3 * M + 1):(4 * M)], M)
  CTSs[[5]][k,] <- make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(4 * M + 1):(5 * M)], M)
  CTSs[[6]][k,] <- make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(5 * M + 1):(6 * M)], M)
  
  
}



#################    ################
# SUMMARIZE ALL THE RESULTS!
#################    ################
hope <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.means)), n = 4))

names(hope) <- c("Bias", "Variance", "Average SD", "95% CI Coverage", "95% CI Length", "root(MSE)")
hope

print(xtable(hope, caption = "Simulation results from K = 2500 iterations", type = "latex"), file = "CTS.summary.tex")


########## Plots for catching zero and nonzero!

spike.summary <- as.data.frame(spike.summary)
spike.summary.nonzero <- as.data.frame(spike.summary.nonzero)

names(spike.summary) <- c("sigma=0.1", "sigma=0.25", "sigma=0.5")
names(spike.summary.nonzero) <- c("sigma=0.1, delta=1", "sigma=0.1, delta=2",
                                  "sigma=0.25, delta=1", "sigma=0.25, delta=2",
                                  "sigma=0.5, delta=1", "sigma=0.5, delta=2")


##### Saving all the simulation iterations 
all.CTS.summary <- as.data.frame(do.call(cbind, CTSs))

save(spike.summary, file = "meta.catch.spike.R")
save(spike.summary.nonzero, file = "meta.catch.nospike.R")
save(all.CTS.summary, file = "meta.CTSs.R")


