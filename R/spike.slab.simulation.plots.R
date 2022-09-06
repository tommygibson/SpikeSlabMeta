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
library(Cairo)
library(latex2exp)


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
dat.box <- cbind.data.frame(c(as.vector(unlist(spike)), as.vector(unlist(spike.summary.nonzero))),
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

spike$sigma.delta <- factor(spike$sigma.delta,
                            levels = c(0.1, 0.25, 0.5),
                            labels = c(TeX("$\\sigma_{\\delta}$ = 0.1"), 
                                       TeX("$\\sigma_{\\delta}$ = 0.25"),
                                       TeX("$\\sigma_{\\delta}$ = 0.5")))

setEPS()
cairo_ps(filename = "TeX/spike_boxplot.eps", height = 4, width = 5)
spike %>%
  ggplot(mapping = aes(x = as.factor(delta0), y = spike)) +
  theme_bw() +
  theme(text = element_text(size = 12, family = "Times New Roman")) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1, alpha = 1) +
  labs(x = TeX("$\\delta_0$"), y = TeX("P($\\delta_0 = 0$ | $Y_k$)")) +
  #labs(x = expression(delta[0]), y = expression(paste('P(', delta[0], '=0 | ', Y[k], ')'))) +
  ylim(c(0, 1)) +
  facet_wrap(~sigma.delta,
             labeller = label_parsed)
dev.off()

# ggsave("TeX/spike_boxplot.pdf", plot = spike_plot,
#        height = 4, width = 5, units = "in")
# ggsave(filename = "TeX/spike_boxplot.eps", plot = spike_plot,
#        height = 4, width = 5, units = "in", 
#        device = cairo_ps, dpi = 600)



