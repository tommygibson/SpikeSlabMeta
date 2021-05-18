########## Show how changing L (# of draws from predictive dist'n in each iteration of MCMC)
########## changes inferences


library(R2jags)
library(here)
library(tidyverse)
library(ggplot2)

source("spike.functions.R")

#### Make plot for a given set of hyperparameters

delta0 <- 2
beta0 <- log(.15 / .85)
nu0 <- log(.15 / .85)

sigma.d <- 0.8
sigma.b <- 0.8
sigma.n <- 0.8

S <- 10
L <- c(1, 10, 100, 1000)

LRp <- matrix(nrow = 10000, ncol = length(L))

set.seed(1515)

delta.i <- rnorm(S, mean = delta0, sd = sigma.d)
beta.i <- rnorm(S, mean = beta0, sd = sigma.b)
nu.i <- rnorm(S, mean = nu0, sd = sigma.n)

pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))

probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities

tables <- t(apply(probs, 1, rmultinom, n = 1, size = N)) # generate contingency tables

y <- tables[,c(1, 3)]
n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
n.tot <- rep(N, S)

for(i in 1:length(L)){
  
  # prepare data for going into jags

  M <- L[i]

  
  meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                   a = -1, b = 0.5, c = 0, d = 0.5, e = -1, f = 0.5)
  
  meta.params <- c("LRpnew")
  
  fit <- jags(data = meta.dat, inits = init.gen, parameters.to.save = meta.params,
              model.file = "R/meta_confusion.txt",
              n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  
  if(M < 2) {
    LRp[, i] <- fit$BUGSoutput$sims.matrix
  }
  else {
    LRp[, i] <- apply(fit$BUGSoutput$sims.matrix, 1, mean)
  }
  
  
}

LRp.dat <- cbind.data.frame(as.vector(LRp),
                            as.factor(rep(L, each = dim(LRp)[1])))

names(LRp.dat) <- c("LR+", "L")

ggplot(data = LRp.dat, mapping = aes(x = `LR+`, linetype = L, size = L)) +
  geom_density(fill = "black", alpha = 0.05) +
  xlim(c(0.5, 8.5)) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_size_manual(values = c(0.5, 0.5, 1, 0.5)) +
  theme_bw() +
  labs(title = "Posterior distribution of LR+ for varying L", ylab = "")

apply(LRp, 2, mean)
apply(LRp, 2, var)
apply(LRp, 2, quantile, c(.025, .975))

ggsave("TeX/L_plot.pdf", plot = last_plot(), device = "pdf",
       width = 5, height = 5, units = "in")




                            