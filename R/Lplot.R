########## Show how changing L (# of draws from predictive dist'n in each iteration of MCMC)
########## changes inferences


library(R2jags)
library(here)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(xtable)

# font_import(pattern = "lmroman*")

source("R/spike.functions.R")

#### Make plot for a given set of hyperparameters

delta0 <- 2
beta0 <- log(.15 / .85)
nu0 <- log(.15 / .85)

sigma.d <- sigma.b <- sigma.n <- c(0.5, 1)

S <- 10
# N <- 500
K <- 100
L <- c(1, 10, 100, 1000, 10000)

LRp <- as.data.frame(matrix(nrow = K * length(L) * length(sigma.d), ncol = 5))
names(LRp) <- c("iteration", "sigma", "L", "SD", "95% Length")

set.seed(1516)

for(j in 1:length(sigma.d)){
  for(k in 1:K){
    
    N <- floor(runif(S, 250, 2501))
    delta.i <- rnorm(S, mean = delta0, sd = sigma.d[j])
    beta.i <- rnorm(S, mean = beta0, sd = sigma.b[j])
    nu.i <- rnorm(S, mean = nu0, sd = sigma.n[j])
    
    pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
    pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
    pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
    pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))
    
    probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities
    
    # generate contingency tables
    tables <- matrix(nrow = S, ncol = 4)
    for(i in 1:S){
      tables[i,] <- rmultinom(1, size = N[i], prob = probs[i,])
    }
    
    y <- tables[,c(1, 3)]
    n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
    n.tot <- N
    
    # prepare data for going into jags
    M <- 1
    B.beta <- B.delta <- B.nu <- 4
    a <- -1
    b <- 0.25
    c <- 0
    d <- 0.25
    e <- -1
    f <- 0.25
    
    meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                     a = a, b = b, c = c, d = d, e = e, f = f,
                     B.beta = B.beta, B.delta = B.delta, B.nu = B.nu)
    
    meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
    
    fit <- do.call(jags.parallel,
                   list(data = names(meta.dat), inits = init.gen, parameters.to.save = meta.params,
                        model.file = "R/Models/meta_confusion.txt",
                        n.chains = 4, n.iter = 3000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
    
    # second column of function output holds LR+
    LRp.1 <- make_cts0_from_hyper(fit$BUGSoutput$sims.matrix, M = 1)[,2]
    LRp.10 <- make_cts0_from_hyper(fit$BUGSoutput$sims.matrix, M = 10)[,2]
    LRp.100 <- make_cts0_from_hyper(fit$BUGSoutput$sims.matrix, M = 100)[,2]
    LRp.1000 <- make_cts0_from_hyper(fit$BUGSoutput$sims.matrix, M = 1000)[,2]
    LRp.10000 <- make_cts0_from_hyper(fit$BUGSoutput$sims.matrix, M = 10000)[,2]
    
    LRp[(K * length(L) * (j - 1)  + (k - 1) * length(L) + 1):(K * length(L) * (j - 1) + (k - 1) * length(L) + length(L)),] <-
      cbind(rep(k, length(L)), rep(sigma.d[j], length(L)), L, rbind(c(sd(LRp.1[LRp.1 < 30]), diff(quantile(LRp.1[LRp.1 < 30], c(.025, .975)))),
                                                                    c(sd(LRp.10[LRp.10 < 30]), diff(quantile(LRp.10[LRp.10 < 30], c(.025, .975)))),
                                                                    c(sd(LRp.100[LRp.100 < 30]), diff(quantile(LRp.100[LRp.100 < 30], c(.025, .975)))),
                                                                    c(sd(LRp.1000[LRp.1000 < 30]), diff(quantile(LRp.1000[LRp.1000 < 30], c(.025, .975)))),
                                                                    c(sd(LRp.10000[LRp.10000 < 30]), diff(quantile(LRp.10000[LRp.10000 < 30], c(.025, .975))))))
  }   
}


L.table <- LRp %>%
  group_by(sigma, L) %>%
  # ggplot(aes(x = `95% Length`, group = interaction(L, sigma))) +
  #   geom_density(aes(color = sigma, linetype = as.factor(L)))
  summarize(mean.sd = mean(SD),
            mean.95 = mean(`95% Length`)) %>%
  pivot_wider(names_from = sigma, values_from = c(mean.sd, mean.95)) %>%
  select(c(1, 2, 4, 3, 5))

print(xtable(L.table, digits = c(0, 0, 3, 3, 3, 3), type = "latex"), file = here("TeX", "L.table.tex"), include.rownames = FALSE)

LRp.low <- cbind.data.frame(as.vector(LRp[[1]]),
                            as.factor(as.character(rep(L, each = dim(LRp[[1]])[1]))))
LRp.high <- cbind.data.frame(as.vector(LRp[[2]]),
                            as.factor(as.character(rep(L, each = dim(LRp[[2]])[1]))))

names(LRp.low) <- names(LRp.high) <- c("LR+", "L")


LRp.lowhigh <- cbind.data.frame(rbind(LRp.low, LRp.high), 
                                as.factor(rep(c("SDs = 0.5", "SDs = 1"), each = dim(LRp.low)[1])))
names(LRp.lowhigh) <- c("LR+", "L", "SDs")

Lplot <- LRp.lowhigh %>%
           ggplot(mapping = aes(x = `LR+`, linetype = L, size = L)) +
             theme_bw() +
             theme(text = element_text(size = 12, family = "LM Roman 10")) +
             geom_density(fill = "black", alpha = 0.05) +
             facet_wrap( ~ SDs) +
             xlim(c(0, 10)) +
             scale_linetype_manual(values = c("dotted", "dotdash", "solid", "dashed")) +
             scale_size_manual(values = c(0.5, 0.5, 0.5, 0.5)) + 
             labs(ylab = c("", ""))


ggsave("TeX/L_plot.pdf", plot = Lplot, device = "pdf",
       width = 5, height = 5, units = "in")

simsum <- function(x){
  sd <- apply(x, 2, sd)
  length <- as.vector(diff(apply(x, 2, quantile, c(.025, .975))))
  
  return(sigfig(cbind(sd, length), n = 3))
}

sd0.5 <- simsum(LRp[[1]])
sd1 <- simsum(LRp[[2]])

Lplot.tab <- cbind.data.frame(sd0.5, sd1)
rownames(Lplot.tab) <- c("L = 1", "L = 10", "L = 100", "L = 1000")

print(xtable(Lplot.tab, caption = "Posterior distribution of LR+ for varying L for two simulated datasets",
             type = "latex"),
      file = "TeX/Lplot.tex")

##################################################
### A different approach
##################################################
### Take a given gamma, see sd(CTS0(gamma)) for varying L
##################################################

## change of plans, don't need different L's
## take L = 1, and the rest are just divided by sqrt(L)
## for a standard error

nsims <- 10000
gamma1 <- c(log(.15 / .85), 2, log(.15 / .85), 0.5, 0.5, 0.5)
gamma2 <- c(log(.15 / .85), 2, log(.15 / .85), 1, 1, 1)

gamma1.mat <- matrix(rep(gamma1, nsims), nrow = nsims, byrow = T)
gamma2.mat <- matrix(rep(gamma2, nsims), nrow = nsims, byrow = T)
gamma.L <- matrix(nrow = nsims, ncol = 2)
#gamma1.L <- gamma2.L <- matrix(nrow = nsims, ncol = 5)

set.seed(613)

gamma.L[,1] <- make_cts0_from_hyper(gamma1.mat, M = 1)[,2]
# gamma1.L[,2] <- make_cts0_from_hyper(gamma1.mat, M = 10)[,2]
# gamma1.L[,3] <- make_cts0_from_hyper(gamma1.mat, M = 100)[,2]                            
# gamma1.L[,4] <- make_cts0_from_hyper(gamma1.mat, M = 1000)[,2]
# gamma1.L[,5] <- make_cts0_from_hyper(gamma1.mat, M = 10000)[,2]

gamma.L[,2] <- make_cts0_from_hyper(gamma2.mat, M = 1)[,2]
# gamma2.L[,2] <- make_cts0_from_hyper(gamma2.mat, M = 10)[,2]
# gamma2.L[,3] <- make_cts0_from_hyper(gamma2.mat, M = 100)[,2]                            
# gamma2.L[,4] <- make_cts0_from_hyper(gamma2.mat, M = 1000)[,2]
# gamma2.L[,5] <- make_cts0_from_hyper(gamma2.mat, M = 10000)[,2]

apply(gamma.L, 2, sd)

# one.gamma.L <- cbind.data.frame(c(1, 10, 100, 1000, 10000),
#                                 apply(gamma1.L, 2, function(x) sd(x[x < 30])),
#                                 apply(gamma2.L, 2, function(x) sd(x[x < 30])))
# names(one.gamma.L) <- c("L", "$sigma_{beta} = sigma_{delta} = sigma_{nu} = 0.5$",
#                         "$sigma_{beta} = sigma_{delta} = sigma_{nu} = 1$")
# 
# print(xtable(one.gamma.L, digits = c(0, 0, 3, 3), type = "latex"), file = here("TeX", "one.gamma.L.tex"), include.rownames = FALSE)
