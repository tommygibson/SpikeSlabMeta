#### taking small results from cts simulation to see how many iterations are needed

library(tidyverse)
library(here)
library(xtable)
library(latex2exp)
library(extrafont)

spike <- readRDS(here("R", "Results", "meta.spike.results.rds"))[[1]]

(spike.table <- spike %>% 
  group_by(delta0, sigma.delta) %>%
  summarize(Mean = mean(spike),
            SD = sd(spike)) %>%
    pivot_longer(!c(delta0, sigma.delta), names_to = "stat", values_to = "count") %>%
    arrange(stat) %>%
    pivot_wider(names_from = delta0, 
                values_from = count) %>%
    relocate(stat))

# print(xtable(spike.table, caption = "Spike table", type = "latex", digits = 4),
#       file = "TeX/spike.table.tex", include.rownames = FALSE)
  
  
spike$sigma.delta <- factor(spike$sigma.delta,
                            levels = c(0.1, 0.25, 0.5),
                            labels = c("sigma[delta] == 0.1", 
                                       "sigma[delta] == 0.25",
                                       "sigma[delta] == 0.5"))

spike.boxplot <- spike %>%
  ggplot(mapping = aes(x = as.factor(delta0), y = spike)) +
  theme_bw() +
  theme(text = element_text(size = 12, family = "LM Roman 10")) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1, alpha = 1) +
  labs(x = TeX("$\\delta_0$"), y = TeX("P($\\delta_0 = 0$ | $Y_k$)")) +
  #labs(x = expression(delta[0]), y = expression(paste('P(', delta[0], '=0 | ', Y[k], ')'))) +
  ylim(c(0, 1)) +
  facet_wrap(~sigma.delta,
             labeller = label_parsed)

# ggsave("TeX", "spike_boxplot.pdf", spike.boxplot,
#        device = "pdf", height = 4, width = 5, units = "in", dpi = 600)
