# simulation to estimate CTSs with 3RE meta-analysis model
library(R2jags)
library(knitr)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(here)

setwd("/u/home/t/tagibson/Projects/SpikeSlabMeta")
i_am("spike.slab.simulation.R")

source(here("spike.functions.R"))
# 
# heart.disease <- read.csv(here('syncope.cleaned.csv')) %>%
#   filter(counts == 1, Variable == "Heart Disease")


####### Simulation: BIAS, COVERAGE %, 95% LENGTH FOR CTSs
# factorial design, S \in (5, 10, 30)
# N will be random each iteration
# S <- dim(heart.disease)[1]
# N <- heart.disease$N_i

S.list <- c(5, 10, 30)

# K determined by MCSE of bias for LR+ 
# Var(\hat{\theta}) \le 0.25, want MCSE < 0.01
set.seed(12) #re-run april 4

# start with 100 iterations to see what variances look like
K <- 1500

sigma.nu <- sigma.beta <- 0.4
nu0 <- beta0 <- log(0.15 / (1 - 0.15))

delta0 <- 2
sigma.delta <- c(.1, .25, .5)

CTSs <- list()
for(i in 1:6){
  CTSs[[i]] <- as.data.frame(matrix(nrow = (2 * K * length(sigma.delta) * length(S.list)), ncol = 9))
  names(CTSs[[i]]) <- c("Method", "Iteration", "S", "sigma.delta", "mean.est", "sd.est", "ci.lower", "median.est", "ci.upper")
}

names(CTSs) <- c("LRm", "LRp", "NPV", "PPV", "sens", "spec")


cts.seeds <- list()
for(i in 3:3){
  for(j in 1:length(sigma.delta)){
    
    for(k in 1:K){
      
      # random study-level sample sizes between 250 and 2500
      cts.seeds[[(i - 1) * 3 * K + (j - 1) * K + k]] <- .Random.seed
      S <- S.list[i]
      N <- floor(runif(S, 250, 2501))
      
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
      
      meta.params <- c("LRmnew", "LRpnew", "sensnew", "specnew", "PPVnew", "NPVnew",
                       "z.LRm", "z.LRp", "z.sens", "z.spec", "z.PPV", "z.NPV")
      if(S <= 5){
        std.model <- do.call(jags.parallel,
                             list(data = names(meta.dat), inits = init.gen.gamma, parameters.to.save = meta.params,
                                  model.file = here("meta_confusion_gamma.txt"),
                                  n.chains = 4, n.iter = 4000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
      } else if(S > 5){
        std.model <- do.call(jags.parallel,
                             list(data = names(meta.dat), inits = init.gen, parameters.to.save = meta.params,
                                  model.file = here("meta_confusion.txt"),
                                  n.chains = 4, n.iter = 4000, n.thin = 2, n.burnin = 2000, DIC = FALSE))
      }
      # monte carlo estimates
      CTSs[[1]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, S, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,1:M], M))
      CTSs[[2]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, S, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(M + 1):(2 * M)], M))
      CTSs[[3]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, S, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(2 * M + 1):(3 * M)], M))
      CTSs[[4]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, S, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(3 * M + 1):(4 * M)], M))
      CTSs[[5]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, S, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(4 * M + 1):(5 * M)], M))
      CTSs[[6]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k - 1)),] <- c("CTS0", k, S, sigma.delta[j], make_CTS_sum1(std.model$BUGSoutput$sims.matrix[,(5 * M + 1):(6 * M)], M))
      # plug-in estimates
      CTSs[[1]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, S, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 1), c(1, 2, 3, 5, 7)])
      CTSs[[2]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, S, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 2), c(1, 2, 3, 5, 7)])
      CTSs[[3]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, S, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 3), c(1, 2, 3, 5, 7)])
      CTSs[[4]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, S, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 4), c(1, 2, 3, 5, 7)])
      CTSs[[5]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, S, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 5), c(1, 2, 3, 5, 7)])
      CTSs[[6]][((i - 1) * 6 * K + (j - 1) * 2 * K + (2 * k)),] <- c("CTSh", k, S, sigma.delta[j], std.model$BUGSoutput$summary[(6 * M + 6), c(1, 2, 3, 5, 7)])
      
      
      
    }
  }
}


# CTSs$LRp %>% 
#   type_convert() %>%
#   group_by(S, sigma.delta, Method) %>%
#   summarize(v = var(mean.est),
#             m = mean(mean.est, na.rm = TRUE),
#             samp = v / .01^2)
cts <- list(CTSs, cts.seeds)
names(cts) <- c("CTSs", "seeds")
saveRDS(cts, here("Results", "cts.sim.30.rds"))
# cts.results <- list(CTSs, cts.seeds)
# names(cts.results) <- c("CTSs", "seeds")
# # hope <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.means)), n = 4))
# # hope.med <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.meds)), n = 4))
# # 
# # names(hope) <- c("Bias", "Variance", "Average SD", "95% CI Coverage", "95% CI Length", "root(MSE)")
# # hope
# # 
# # print(xtable(hope, caption = "Simulation results from K = 2500 iterations", type = "latex"), file = "TeX/CTS.summary.tex")
# 
# saveRDS(cts.results, here("Results", "cts.simulation.results.rds"))
