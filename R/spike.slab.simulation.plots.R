############ MAKING PLOTS AND TABLES OUT OF SPIKE/SLAB SIMULATION RESULTS

setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta/Code")

library(knitr)
library(gplots)
library(pander)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(gridExtra)

font_import(pattern = "lmroman*")

##### LOAD DATA FROM SIMULATIONS
load("meta.catch.spike.R")
load("meta.catch.nospike.R")
load("meta.CTSs.R")


########### FUNCITONS WE'LL USE
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
  sd.est <- sd(x[,1]) # mean of estimator
  avg.sd <- mean(x[,2]) # average SD
  cover.95 <- sum(x[,3] < target[1] & x[,5] > target[1]) / dim(x)[1] # 95% coverage probability
  length.95 <- mean(x[, 5] - x[, 3]) # average 95% length
  RMSE <- sqrt(mean((x[, 1] - target[1])^2)) # root(MSE)
  MCSE <- sqrt(var(x[,1]) / dim(x)[1])  # Monte Carlo standard error
  
  return(c(bias, sd.est, avg.sd, cover.95, length.95, RMSE, MCSE))
}

reject_by_cutoff <- function(cut, mat){
  
  return(apply(mat, 2, function(x) sum(x > cut) / length(x)))
  
}



#### make summaries to actually put in the paper

############# SIMULATIONS 1 AND 2: SPIKE/SLAB

# for table summary and boxplot
dat.box <- cbind.data.frame(c(as.vector(unlist(spike.summary)), as.vector(unlist(spike.summary.nonzero))),
                            c(rep(c("SD = 0.1", "SD = 0.25", "SD = 0.5"), each = 1000), rep(c("SD = 0.1", "SD = 0.25", "SD = 0.5"), each = 2000)),
                            c(rep("0", 3000), rep(rep(c("1", "2"), each = 1000), 3)))
names(dat.box) <- c("P(=0)", "sigma", "delta0")

spike.table.data <- as.data.frame(matrix(aggregate(dat.box$`P(=0)` ~ sigma + delta0, FUN = mean, data = dat.box)[,3], 
                                         nrow = 3, byrow = FALSE))
print(xtable(spike.table.data, caption = "Average spike height for each simulation scenario", type = "latex", digits = c(0, 4, 4, 4)), 
      file = "spike.summary.tex", include.rownames = TRUE)


######## boxplots of zero/nonzero spikes

spike_plot <- dat.box %>%
  ggplot(mapping = aes(x = delta0, y = `P(=0)`)) +
  theme_bw() +
  theme(text = element_text(size = 12, family = "LM Roman 10")) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1, alpha = 1) +
  labs(x = "Mean of REs", y = "Height of spike") +
  ylim(c(0, 1)) +
  facet_wrap(~sigma)

ggsave("spike_boxplot.pdf", plot = spike_plot, 
       path = "/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta",
       height = 4, width = 5, units = "in")

####### Plots for catching true zeros and nonzeros

catch.zero <- cbind.data.frame(as.vector(t(sapply(seq(0.01, 1, 0.01), reject_by_cutoff, mat = spike.summary))),
                               seq(0.01, 1, 0.01),
                               rep(names(spike.summary), each = 100))
names(catch.zero) <- c("P.rej", "cutoff", "situation")


zero.plot <- ggplot(catch.zero, mapping = aes(x = cutoff, y = P.rej, linetype = situation)) + 
  geom_line(size = 0.75) +
  theme_bw() +
  theme(text = element_text(size = 10, family = "LM Roman 10")) +
  labs(title = "Zero mean effect",
       x = "Cutoff value", y = "P(Spike > Cutoff)", linetype = "SD of REs") +
  scale_linetype_manual(values = c(1, 2, 3)) + 
  #geom_hline(yintercept = 0.95, color = "gray65", linetype = 6) +
  scale_y_continuous(breaks = seq(0, 1, 0.2))

catch.nonzero <- cbind.data.frame(as.vector(t(sapply(seq(0.01, 1, 0.01), reject_by_cutoff, 
                                                     mat = spike.summary.nonzero))),
                                  seq(0.01, 1, 0.01),
                                  rep(names(spike.summary.nonzero), each = 100),
                                  c(rep("0.1", 200), rep("0.25", 200), rep("0.5", 200)),
                                  rep(rep(c("1", "2"), each = 100), 3))

names(catch.nonzero) <- c("P.accept", "cutoff", "situation", "sigma", "delta0")

nonzero.plot <- ggplot(catch.nonzero, mapping = aes(x = cutoff, y = P.accept, linetype = sigma, size = delta0)) + 
  geom_line() +
  theme_bw() + 
  theme(text = element_text(size = 10, family = "LM Roman 10")) +
  labs(title = "Non-zero mean effect",
       x = "Cutoff value", y = "P(Spike > Cutoff)", size = "Mean of REs", linetype = "SD of REs") +
  scale_size_manual(values = c(0.5, 1)) + 
  scale_linetype_manual(values = c(1, 2, 3)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) 

SSP.plot <- grid.arrange(zero.plot, nonzero.plot, nrow = 2)

ggsave("SSP-plot.pdf", plot = SSP.plot, 
       path = "/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta",
       width = 6, height = 6, units = "in")
# ggsave("zero-plot.pdf", plot = zero.plot,
#        path = "/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta",
#        width = 6, height = 3, units = "in")


######### SIMULATION 3: CTSs
CTSs <- list()

for(i in 1:6){
  CTSs[[i]] <- all.CTS.summary[,(5 * (i - 1) + 1):(5 * i)]
}

#### Expected values of parameters

# FIND EXPECTED VALUES FOR STATS GIVEN HYPERPARAMETER VALUES

nsim <- 50000
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


hope <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.means)), n = 4))

names(hope) <- c("Bias", "Variance", "Average SD", "95% CI Coverage", "95% CI Length", "root(MSE)", "MCSE(bias)")
hope

print(xtable(hope, caption = "Results from Simulation 3 with K_3 = 2500 repetitions", type = "latex"), file = "CTS.summary.04.19.2021.tex")
