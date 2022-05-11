#### testing post-processing speed vs within-MCMC speed
### do we get the same results by just following hyperparameters and estimating CTS0 afterwards?

set.seed(1234)

library(R2jags)
library(knitr)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(here)

source(here("R", "spike.functions.R"))

S <- 5
N <- floor(runif(10, 250, 2501))
sigma.nu <- sigma.beta <- 0.4
nu0 <- beta0 <- log(0.15 / (1 - 0.15))
delta0 <- 2
delta0fix <- 0
sigma.delta <- 0.5

delta0 <- 2

#### generate a single dataset to see if we can recreate CTS0 from
#### just the hyperparameters and post-processing

delta.i <- rnorm(S, mean = delta0, sd = sigma.delta)
beta.i <- rnorm(S, mean = beta0, sd = sigma.beta)
nu.i <- rnorm(S, mean = nu0, sd = sigma.nu)

pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))

probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities

tables <- matrix(0, nrow = S, ncol = 4)
for(s in 1:S){
 tables[s,] <- rmultinom(probs[s,], n = 1, size = N[i]) # generate contingency tables
}

# prepare data for going into jags
y <- tables[,c(1, 3)]
n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
n.tot <- N
M <- 100

a <- -1
b <- 0.5
c <- 0
d <- 0.5
e <- -1
f <- 0.5

meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                a = a, b = b, c = c, d = d, e = e, f = f)

params.hyper <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
params.cts <- c("LRmnew", "LRpnew", "sensnew", "specnew", "PPVnew", "NPVnew",
                "z.LRm", "z.LRp", "z.sens", "z.spec", "z.PPV", "z.NPV")

fit.cts <- do.call(jags.parallel,
                    list(data = names(meta.dat), inits = init.gen, parameters.to.save = params.cts,
                         model.file = here("R", "meta_confusion.txt"),
                         n.chains = 4, n.iter = 4000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
fit.hyper <- do.call(jags.parallel,
                     list(data = names(meta.dat), inits = init.gen, parameters.to.save = params.hyper,
                          model.file = here("R", "meta_confusion.txt"),
                          n.chains = 4, n.iter = 4000, n.thin = 2, n.burnin = 2000, DIC = FALSE))

LRm_0 <- make_CTS_sum1(fit.cts$BUGSoutput$sims.matrix[,1:M], M)
LRp_0 <- make_CTS_sum1(fit.cts$BUGSoutput$sims.matrix[,(M + 1):(2 * M)], M)
NPV_0 <- make_CTS_sum1(fit.cts$BUGSoutput$sims.matrix[,(2 * M + 1):(3 * M)], M)
PPV_0 <- make_CTS_sum1(fit.cts$BUGSoutput$sims.matrix[,(3 * M + 1):(4 * M)], M)
sens_0 <- make_CTS_sum1(fit.cts$BUGSoutput$sims.matrix[,(4 * M + 1):(5 * M)], M)
spec_0 <- make_CTS_sum1(fit.cts$BUGSoutput$sims.matrix[,(5 * M + 1):(6 * M)], M)

LRm_h <- fit.cts$BUGSoutput$summary[(6 * M + 1), c(1, 2, 3, 5, 7)]
LRp_h <- fit.cts$BUGSoutput$summary[(6 * M + 2), c(1, 2, 3, 5, 7)]
NPV_h <- fit.cts$BUGSoutput$summary[(6 * M + 3), c(1, 2, 3, 5, 7)]
PPV_h <- fit.cts$BUGSoutput$summary[(6 * M + 4), c(1, 2, 3, 5, 7)]
sens_h <- fit.cts$BUGSoutput$summary[(6 * M + 5), c(1, 2, 3, 5, 7)]
spec_h <- fit.cts$BUGSoutput$summary[(6 * M + 6), c(1, 2, 3, 5, 7)]

hyper.full <- fit.hyper$BUGSoutput$sims.matrix


apply(make_cts0_from_hyper(hyper.full, M = 100), 2, mean)
t(matrix(c(LRm_0, LRp_0, NPV_0, PPV_0, sens_0, spec_0), nrow = 6, byrow = T)[,1])
apply(make_ctsh_from_hyper(hyper.full), 2, mean)
t(matrix(c(LRm_h, LRp_h, NPV_h, PPV_h, sens_h, spec_h), nrow = 6, byrow = T)[,1])


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
    