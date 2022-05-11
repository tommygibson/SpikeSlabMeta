######### simulation for spike/slab meta-analysis

set.seed(1234)

library(R2jags)
library(knitr)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(here)

source(here("R", "spike.functions.R"))

heart.disease <- read.csv(here('R', 'syncope.cleaned.csv')) %>%
  filter(counts == 1, Variable == "Heart Disease")

# next line only necessary if hasn't been imported before
# font_import(pattern = "lmroman*")



########### DATA GENERATION

# FIND EXPECTED VALUES FOR STATS GIVEN HYPERPARAMETER VALUES

nsim <- 10000
sigma.nu <- sigma.beta <- 0.4
nu0 <- beta0 <- log(0.15 / (1 - 0.15))
delta0 <- 2
delta0fix <- 0
sigma.delta <- 0.5

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
stats.meds <- apply(stats.sims, 2, quantile, 0.5)


###### Okay, now actual simulation
###### need to generate 2x2 tables and then apply models to the data repeatedly


###### Spike/slab model
###### Assess how well it catches true zeros and avoids real nonzeros

#### FIRST, HOW WELL IT CATCHES ZEROS

# get number of studies and sample size per study from syncope heart disease
S <- dim(heart.disease)[1]
N <- heart.disease$N_i
# K <- 1000
K <- 5
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
    
    tables <- matrix(0, nrow = S, ncol = 4)
    for(j in 1:S){
      tables[j,] <- rmultinom(probs[j,], n = 1, size = N[j]) # generate contingency tables
    }
    
    # prepare data for going into jags
    y <- tables[,c(1, 3)]
    n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
    n.tot <- rep(N, S)
    
    M <- 100
    a <- -1
    b <- 0.5
    c <- 0
    d <- 0.5
    e <- -1.5
    f <- 0.5
    p <- 0.5
    
    meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                     a = a, b = b, c = c, d = d, e = e, f = f, p = p)
    
    
    # as a first run we'll just follow the hyperparameters
    meta.params <- c("spike")
    
    spike.zero <- do.call(jags.parallel,
                          list(data = names(meta.dat), inits = init.gen.spike, parameters.to.save = meta.params,
                          model.file = "R/meta_confusion_spike.txt",
                          n.chains = 4, n.iter = 5000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
    
    # columns 1 and 3 are mean (1 - rho)
    # columns 2 and 4 are indicators for mean (1 - rho) < x
    spike.summary[k, i] <- spike.zero$BUGSoutput$summary[1]
    
    
  }
}

seed.after.zeros <- .Random.seed

####### SECOND: HOW WELL IT CATCHES NON-ZEROS


#K <- 1000
K <- 5
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
      
      tables <- matrix(0, nrow = S, ncol = 4)
      for(l in 1:S){
        tables[l,] <- rmultinom(probs[l,], n = 1, size = N[l]) # generate contingency tables
      }
      
      # prepare data for going into jags
      y <- tables[,c(1, 3)]
      n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
      n.tot <- rep(N, S)
      M <- 100
      a <- -1
      b <- 0.5
      c <- 0
      d <- 0.5
      e <- -1.5
      f <- 0.5
      p <- 0.5
      
      meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                       a = a, b = b, c = c, d = d, e = e, f = f, p = p)
      
      # as a first run we'll just follow the hyperparameters
      meta.params <- c("spike")
      
      spike.zero <- do.call(jags.parallel,
                            list(data = names(meta.dat), inits = init.gen.spike, parameters.to.save = meta.params,
                            model.file = "R/meta_confusion_spike.txt",
                            n.chains = 4, n.iter = 5000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
      
      # columns represent the two values for sigma_delta (small and moderate)
      
      spike.summary.nonzero[k, (length(delta0.nonzero) * (i - 1) + j)] <- spike.zero$BUGSoutput$summary[1]
      
      
    }
  }
}

seed.after.nonzero <- .Random.seed



####### THIRD: BIAS, COVERAGE %, 95% LENGTH FOR CTSs

# K determined by MCSE of bias for LR+ 
# Var(\hat{\theta}) \le 0.25, want MCSE < 0.01
#K <- 2500
K <- 5

delta0 <- 2
sigma.delta <- c(.1, .25, .5)

CTSs <- list()
for(i in 1:6){
  CTSs[[i]] <- as.data.frame(matrix(nrow = (2 * K * length(sigma.delta)), ncol = 8))
  names(CTSs[[i]]) <- c("Method", "Iteration", "sigma.delta", "mean.est", "sd.est", "ci.lower", "median.est", "ci.upper")
}

names(CTSs) <- c("LRm", "LRp", "NPV", "PPV", "sens", "spec")


# .Random.seed <- seed.after.nonzero
for(j in 1:length(sigma.delta)){
  
  for(k in 1:K){
      
    delta.i <- rnorm(S, mean = delta0, sd = sigma.delta[j])
    beta.i <- rnorm(S, mean = beta0, sd = sigma.beta)
    nu.i <- rnorm(S, mean = nu0, sd = sigma.nu)
    
    pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
    pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
    pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
    pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))
    
    probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities
    
    tables <- matrix(0, nrow = S, ncol = 4)
    for(s in 1:S){
      tables[s,] <- rmultinom(probs[s,], n = 1, size = N[s]) # generate contingency tables
    }
    
    # prepare data for going into jags
    y <- tables[,c(1, 3)]
    n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
    n.tot <- rep(N, S)
    M <- 100
    
    a <- -1
    b <- 0.5
    c <- 0
    d <- 0.5
    e <- -1
    f <- 0.5
    
    meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                     a = a, b = b, c = c, d = d, e = e, f = f)
  
    meta.params <- c("LRmnew", "LRpnew", "sensnew", "specnew", "PPVnew", "NPVnew",
                     "z.LRm", "z.LRp", "z.sens", "z.spec", "z.PPV", "z.NPV")
    
    std.model <- do.call(jags.parallel,
                         list(data = names(meta.dat), inits = init.gen, parameters.to.save = meta.params,
                         model.file = here("R", "meta_confusion.txt"),
                         n.chains = 4, n.iter = 4000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
    
    # monte carlo estimates
    CTSs[[1]][((j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,1:M], M))
    CTSs[[2]][((j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(M + 1):(2 * M)], M))
    CTSs[[3]][((j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(2 * M + 1):(3 * M)], M))
    CTSs[[4]][((j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(3 * M + 1):(4 * M)], M))
    CTSs[[5]][((j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(4 * M + 1):(5 * M)], M))
    CTSs[[6]][((j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(5 * M + 1):(6 * M)], M))
    # plug-in estimates
    CTSs[[1]][((j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 1), c(1, 2, 3, 5, 7)])
    CTSs[[2]][((j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 2), c(1, 2, 3, 5, 7)])
    CTSs[[3]][((j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 3), c(1, 2, 3, 5, 7)])
    CTSs[[4]][((j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 4), c(1, 2, 3, 5, 7)])
    CTSs[[5]][((j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 5), c(1, 2, 3, 5, 7)])
    CTSs[[6]][((j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 6), c(1, 2, 3, 5, 7)])
    
    
    
  }
}



#################    ################
# SUMMARIZE ALL THE RESULTS!
#################    ################
hope <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.means)), n = 4))
hope.med <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.meds)), n = 4))

saveRDS(list(hope))
names(hope) <- c("Bias", "Variance", "Average SD", "95% CI Coverage", "95% CI Length", "root(MSE)")
hope

print(xtable(hope, caption = "Simulation results from K = 2500 iterations", type = "latex"), file = "TeX/CTS.summary.tex")


########## Plots for catching zero and nonzero!

spike.summary <- as.data.frame(spike.summary)
spike.summary.nonzero <- as.data.frame(spike.summary.nonzero)

names(spike.summary) <- c("sigma=0.1", "sigma=0.25", "sigma=0.5")
names(spike.summary.nonzero) <- c("sigma=0.1, delta=1", "sigma=0.1, delta=2",
                                  "sigma=0.25, delta=1", "sigma=0.25, delta=2",
                                  "sigma=0.5, delta=1", "sigma=0.5, delta=2")


##### Saving all the simulation iterations 
# all.CTS.summary <- as.data.frame(do.call(cbind, CTSs))

save(spike.summary, file = "R/meta.catch.spike.8.25.R")
save(spike.summary.nonzero, file = "R/meta.catch.nospike.8.25.R")
saveRDS(CTSs, file = "R/meta.CTSs.8.25.R")


