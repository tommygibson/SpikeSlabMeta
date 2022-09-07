# FIND EXPECTED VALUES FOR STATS GIVEN HYPERPARAMETER VALUES

library(here)
# setwd("/u/home/t/tagibson/Projects/SpikeSlabMeta")
# i_am("calculate.targets.R")

nsim <- 500000
sigma.nu <- sigma.beta <- 0.4
nu0 <- beta0 <- log(0.15 / (1 - 0.15))
delta0 <- 2
sigma.delta <- c(.1, .25, .5)
targets <- as.data.frame(matrix(nrow = 6 * length(sigma.delta), ncol = 4))
names(targets) <- c("CTS", "sigma.delta", "target.mean", "target.median")
targets[,1] <- rep(c("LRm", "LRp", "NPV", "PPV", "Sens", "Spec"), length(sigma.delta))
targets[,2] <- rep(sigma.delta, each = 6)

pi1.h <- 1 / (1 + exp(-(beta0 + delta0 / 2)))
pi0.h <- 1 / (1 + exp(-(beta0 - delta0 / 2)))
psi.h <- 1 / (1 + exp(-nu0))

pi11.h <- pi1.h * psi.h
pi10.h <- (1 - pi1.h) * psi.h
pi01.h <- pi0.h * (1 - psi.h)
pi00.h <- (1 - pi0.h) * (1 - psi.h)

sens.h <- pi11.h / (pi11.h + pi01.h)
spec.h <- pi00.h / (pi00.h + pi10.h)

LRp.h <- sens.h / (1 - spec.h)
LRm.h <- (1 - sens.h) / spec.h

for(i in 1:length(sigma.delta)){
  
  deltasims <- rnorm(nsim, mean = delta0, sd = sigma.delta[i])
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
  
  # mean and median targets
  targets[(i * 6 - 5):(i * 6), 3] <- stats.means
  targets[(i * 6 - 5):(i * 6), 4] <- stats.meds
  
}

saveRDS(targets, here("R", "targets.rds"))

