############ MAKING PLOTS AND TABLES OUT OF SPIKE/SLAB SIMULATION RESULTS

library(knitr)
library(tidyverse)
library(xtable)
library(pROC)
library(extrafont)
library(gridExtra)
library(here)
library(ggplot2)
library(scales)
library(here)

library(extrafont)

# remotes::install_version("Rttf2pt1", version = "1.3.8")
# extrafont::font_import()

# This only necessary if it hasn't been done before
# font_import(pattern = "lmroman*")


##### LOAD DATA FROM SIMULATIONS

spike <- readRDS(here("R", "Results", "meta.spike.results.rds"))[[1]]
source(here("R", "spike.functions.R"))

############# Spike/slab simulation

# for table summary and boxplot
dat.box <- cbind.data.frame(c(as.vector(unlist(spike.summary)), as.vector(unlist(spike.summary.nonzero))),
                            c(rep(c("sigma[delta] == 0.1", "sigma[delta] == 0.25", "sigma[delta] == 0.5"), each = 1000), rep(c("sigma[delta] == 0.1", "sigma[delta] == 0.25", "sigma[delta] == 0.5"), each = 2000)),
                            c(rep("0", 3000), rep(rep(c("1", "2"), each = 1000), 3)))
names(dat.box) <- c("P(=0)", "sigma", "delta0")

spike.table.data <- as.data.frame(matrix(aggregate(dat.box$`P(=0)` ~ sigma + delta0, FUN = mean, data = dat.box)[,3], 
                                         nrow = 3, byrow = FALSE))
print(xtable(spike.table.data, caption = "Average spike height for each simulation scenario", type = "latex", digits = c(0, 4, 4, 4)), 
      file = "spike.summary.tex", include.rownames = TRUE)


######## boxplots of zero/nonzero spikes

sigma.labeller <- c(expression(sigma[delta] == 0.1), expression(sigma[delta] == 0.25), expression(sigma[delta] == 0.5))
names(sigma.labeller) <- c(0.1, 0.25, 0.5)

(spike_plot <- spike %>%
  ggplot(mapping = aes(x = as.factor(delta0), y = spike)) +
  theme_bw() +
  theme(text = element_text(size = 12, family = "LM Roman 10")) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1, alpha = 1) +
  labs(x = expression(delta[0]), y = expression(paste('P(', delta[0], '=0 | ', Y[k], ')'))) +
  ylim(c(0, 1)) +
  facet_wrap(~as.factor(sigma.delta)))

ggsave("TeX/spike_boxplot.pdf", plot = spike_plot,
       height = 4, width = 5, units = "in")

####### Plots for catching true zeros and nonzeros

catch.zero <- cbind.data.frame(as.vector(t(sapply(seq(0.01, 1, 0.01), reject_by_cutoff, mat = spike.summary))),
                               seq(0.01, 1, 0.01),
                               rep(c("0.10", "0.25", "0.50"), each = 100))
names(catch.zero) <- c("P.rej", "cutoff", "situation")


zero.plot <- ggplot(catch.zero, mapping = aes(x = cutoff, y = P.rej, linetype = situation)) + 
  geom_line(size = 0.75) +
  theme_bw() +
  theme(text = element_text(size = 10, family = "LM Roman 10")) +
  labs(title = "Zero mean effect",
       x = "Cutoff", 
       y = expression(paste('Proportion of P(', delta[0], '=0 | ', Y[k], ') > Cutoff')), 
       linetype = expression(sigma[delta])) +
  scale_linetype_manual(values = c(1, 2, 3)) + 
  #geom_hline(yintercept = 0.95, color = "gray65", linetype = 6) +
  scale_y_continuous(breaks = seq(0, 1, 0.2))

catch.nonzero <- cbind.data.frame(as.vector(t(sapply(seq(0.01, 1, 0.01), reject_by_cutoff, 
                                                     mat = spike.summary.nonzero))),
                                  seq(0.01, 1, 0.01),
                                  rep(names(spike.summary.nonzero), each = 100),
                                  c(rep("0.10", 200), rep("0.25", 200), rep("0.50", 200)),
                                  rep(rep(c("1", "2"), each = 100), 3))

names(catch.nonzero) <- c("P.accept", "cutoff", "situation", "sigma", "delta0")

nonzero.plot <- ggplot(catch.nonzero, mapping = aes(x = cutoff, y = P.accept, linetype = sigma, size = delta0)) + 
  geom_line() +
  theme_bw() + 
  theme(text = element_text(size = 10, family = "LM Roman 10")) +
  labs(title = "Non-zero mean effect",
       x = "Cutoff", 
       y = expression(paste('Proportion of P(', delta[0], '=0 | ', Y[k], ') > Cutoff')), 
       size = expression(delta[0]), 
       linetype = expression(sigma[delta])) +
  scale_size_manual(values = c(0.5, 1)) + 
  scale_linetype_manual(values = c(1, 2, 3)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) 

SSP.plot <- grid.arrange(zero.plot, nonzero.plot, nrow = 2)

ggsave("TeX/SSP-plot.pdf", plot = SSP.plot,
       width = 6, height = 6, units = "in")
# ggsave("zero-plot.pdf", plot = zero.plot,
#        path = "/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta",
#        width = 6, height = 3, units = "in")


######### SIMULATION 3: CTSs
# CTSs <- list()
# 
# for(i in 1:6){
#   CTSs[[i]] <- all.CTS.summary[,(5 * (i - 1) + 1):(5 * i)]
# }

#### Expected values of parameters

# # FIND EXPECTED VALUES FOR STATS GIVEN HYPERPARAMETER VALUES
# 
# nsim <- 50000
# sigma.nu <- sigma.beta <- 0.4
# nu0 <- beta0 <- log(0.15 / (1 - 0.15))
# delta0 <- 2
# delta0fix <- 0
# sigma.delta <- 0.4
# 
# deltasims <- rnorm(nsim, mean = delta0, sd = sigma.delta)
# nusims <- rnorm(nsim, mean = nu0, sd = sigma.nu)
# betasims <- rnorm(nsim, mean = beta0, sd = sigma.beta)
# 
# 
# pi11sims <- 1 / (1 + exp(-(betasims + deltasims / 2))) * 1 / (1 + exp(-nusims))
# pi10sims <- (1 - 1 / (1 + exp(-(betasims + deltasims / 2)))) * 1 / (1 + exp(-nusims))
# pi01sims <- 1 / (1 + exp(-(betasims - deltasims / 2))) * (1 - 1 / (1 + exp(-nusims)))
# pi00sims <- (1 - 1 / (1 + exp(-(betasims - deltasims / 2)))) * (1 - 1 / (1 + exp(-nusims)))
# 
# pi11fix <- 1 / (1 + exp(-(betasims + delta0fix / 2))) * 1 / (1 + exp(-nusims))
# pi10fix <- (1 - 1 / (1 + exp(-(betasims + delta0fix / 2)))) * 1 / (1 + exp(-nusims))
# pi01fix <- 1 / (1 + exp(-(betasims - delta0fix / 2))) * (1 - 1 / (1 + exp(-nusims)))
# pi00fix <- (1 - 1 / (1 + exp(-(betasims - delta0fix / 2)))) * (1 - 1 / (1 + exp(-nusims)))
# 
# senssims <- pi11sims / (pi11sims + pi01sims)
# specsims <- pi00sims / (pi00sims + pi10sims)
# 
# LRpsims <- senssims / (1 - specsims)
# LRmsims <- (1 - senssims) / specsims
# 
# PPVsims <- pi11sims / (pi11sims + pi10sims)
# NPVsims <- pi00sims / (pi00sims + pi01sims)
# 
# sensfix <- pi11fix / (pi11fix + pi01fix)
# specfix <- pi00fix / (pi00fix + pi10fix)
# 
# LRpfix <- sensfix / (1 - specfix)
# LRmfix <- (1 - sensfix) / specfix
# 
# PPVfix <- pi11fix / (pi11fix + pi10fix)
# NPVfix <- pi00fix / (pi00fix + pi01fix)
# 
# RDfix <- PPVfix - (1 - NPVfix)
# RRfix <- PPVfix / (1 - NPVfix)
# 
# stats.sims <- cbind.data.frame(LRmsims, LRpsims, NPVsims, PPVsims, senssims, specsims)
# stats.fix <- cbind.data.frame(LRmfix, LRpfix, NPVfix, PPVfix, sensfix, specfix)
# names(stats.sims) <- c("LRm", "LRp", "NPV", "PPV", "Sens", "Spec")
# names(stats.fix) <- c("LRm", "LRp", "NPV", "PPV", "Sens", "Spec")
# 
# apply(stats.sims, 2, mean)
# apply(stats.sims, 2, quantile, c(.025, .5, .975))
# 
# apply(stats.fix, 2, mean)
# apply(stats.fix, 2, quantile, c(.025, .5, .975))
# 
# stats.means <- apply(stats.sims, 2, mean)
# 
# 
# hope <- as.data.frame(sigfig(t(mapply(CTS.overall.sum, CTSs, stats.means)), n = 4))
# 
# names(hope) <- c("Bias", "Variance", "Average SD", "95% CI Coverage", "95% CI Length", "root(MSE)", "MCSE(bias)")
# hope
# 
# print(xtable(hope, caption = "Results from Simulation 3 with K_3 = 2500 repetitions", type = "latex"), file = "CTS.summary.04.19.2021.tex")
# 
# 
# ####### Results from redone simulation april 6 2022
# 
# cts <- lapply(cts, function(x) type_convert(x, guess_integer = TRUE))
# 
# cts$LRm <- left_join(cts$LRm, targets %>% filter(CTS == "LRm"), by = "sigma.delta") %>% type_convert(guess_integer = TRUE)
# cts$LRp <- left_join(cts$LRp, targets %>% filter(CTS == "LRp"), by = "sigma.delta") %>% type_convert(guess_integer = TRUE)
# cts$NPV <- left_join(cts$NPV, targets %>% filter(CTS == "NPV"), by = "sigma.delta") %>% type_convert(guess_integer = TRUE)
# cts$PPV <- left_join(cts$PPV, targets %>% filter(CTS == "PPV"), by = "sigma.delta") %>% type_convert(guess_integer = TRUE)
# cts$sens <- left_join(cts$sens, targets %>% filter(CTS == "Sens"), by = "sigma.delta") %>% type_convert(guess_integer = TRUE)
# cts$spec <- left_join(cts$spec, targets %>% filter(CTS == "Spec"), by = "sigma.delta") %>% type_convert(guess_integer = TRUE)
# 
# LRm.summary <- cts$LRm %>%
#   group_by(Method, sigma.delta) %>%
#   summarize(bias = mean(mean.est - target.mean))
# LRp.summary <- cts$LRp %>%
#   group_by(Method, sigma.delta) %>%
#   summarize(bias = mean(mean.est - target.mean),
#             coverage = mean(ci.lower < target.mean & ci.upper > target.mean),
#             rmse = )
# NPV.summary <- cts$NPV %>%
#   group_by(Method, sigma.delta) %>%
#   summarize(bias = mean(mean.est - target.mean))
# PPV.summary <- cts$PPV %>%
#   group_by(Method, sigma.delta) %>%
#   summarize(bias = mean(mean.est - target.mean))
# sens.summary <- cts$sens %>%
#   group_by(Method, sigma.delta) %>%
#   summarize(bias = mean(mean.est - target.mean))
# spec.summary <- cts$spec %>%
#   group_by(Method, sigma.delta) %>%
#   summarize(bias = mean(mean.est - target.mean))



