###### Smaller spike/slab simulation
###### comparing plug-in estimator to monte carlo estimator for one LARGE dataset

set.seed(2021)

library(R2jags)
library(knitr)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(here)

source(here("R", "spike.functions.R"))

S <- 5000

N <- 5000

# means
beta0 <- log(.15 / (1 - .15))
delta0 <- 2
nu0 <- log(.15 / (1 - .15))
# standard deviations
sigma.beta <- 0.75
sigma.delta <- 0.75
sigma.nu <- 0.75

betai <- rnorm(S, beta0, sigma.beta)
deltai <- rnorm(S, delta0, sigma.delta)
nui <- rnorm(S, nu0, sigma.nu)

pi1 <- 1 / (1 + exp(-(betai + deltai / 2)))
pi0 <- 1 / (1 + exp(-(betai - deltai/2)))
psi <- 1 / (1 + exp(-nui))

pi11.i <- pi1 * psi
pi10.i <- (1 - pi1) * psi
pi01.i <- pi0 * (1 - psi)
pi00.i <- (1 - pi0) * (1 - psi)

### for calculating expected statistics vs actual

senssims <- pi11.i / (pi11.i + pi01.i)
specsims <- pi00.i / (pi00.i + pi10.i)

LRpsims <- senssims / (1 - specsims)
LRmsims <- (1 - senssims) / specsims

PPVsims <- pi11.i / (pi11.i + pi10.i)
NPVsims <- pi00.i / (pi00.i + pi01.i)

stats.sims <- cbind.data.frame(LRmsims, LRpsims, NPVsims, PPVsims, senssims, specsims)
sims.means <- apply(stats.sims, 2, mean)

probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities

tables <- t(apply(probs, 1, rmultinom, n = 1, size = N)) # generate contingency tables


# prepare data for going into jags
y <- tables[,c(1, 3)]
n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
n.tot <- rep(N, S)


###### Calculating ACTUAL STATS FOR EACH TABLE

sens.i <- tables[,1] / (tables[,1] + tables[,3])
spec.i <- tables[,4] / (tables[,4] + tables[,2])

LRp.i <- sens.i / (1 - spec.i)
LRm.i <- (1 - sens.i) / spec.i

PPV.i <- tables[,1] / (tables[,1] + tables[,2])
NPV.i <- tables[,4] / (tables[,4] + tables[,3])

stats.real <- cbind.data.frame(LRm.i, LRp.i, NPV.i, PPV.i, sens.i, spec.i)
real.means <- apply(stats.real, 2, mean)

M <- 100

meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                 a = -1, b = 0.5, c = 0, d = 0.5, e = -1, f = 0.5)

meta.params <- c("LRmnew", "LRpnew", "sensnew", "specnew", "PPVnew", "NPVnew", "z.LRm", "z.LRp", "z.NPV", "z.PPV", "z.sens", "z.spec")

std.model <- jags(data = meta.dat, inits = init.gen, parameters.to.save = meta.params,
                  model.file = here("R", "meta_confusion.txt"),
                  n.chains = 4, n.iter = 6000, n.thin = 2, n.burnin = 1000, DIC = FALSE)

mcmc.sims <- std.model$BUGSoutput$sims.matrix
montecarlo.est <- matrix(nrow = 10000, ncol = 6)
montecarlo.est[,1] <- apply(mcmc.sims[,1:M], 1, mean)
montecarlo.est[,2] <- apply(mcmc.sims[,(M + 1):(2 * M)], 1, mean)
montecarlo.est[,3] <- apply(mcmc.sims[,(2 * M + 1):(3 * M)], 1, mean)
montecarlo.est[,4] <- apply(mcmc.sims[,(3 * M + 1):(4 * M)], 1, mean)
montecarlo.est[,5] <- apply(mcmc.sims[,(4 * M + 1):(5 * M)], 1, mean)
montecarlo.est[,6] <- apply(mcmc.sims[,(5 * M + 1):(6 * M)], 1, mean)

plugin.est <- mcmc.sims[,601:606]

apply(montecarlo.est, 2, mean)
apply(plugin.est, 2, mean)





