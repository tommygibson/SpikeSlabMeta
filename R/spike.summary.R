#### taking small results from cts simulation to see how many iterations are needed

library(tidyverse)
library(here)
library(xtable)

spike <- readRDS(here("R", "Results", "meta.spike.results.rds"))[[1]]

(spike.table <- spike %>% 
  filter(delta0 != 0.5) %>%
  group_by(delta0, sigma.delta) %>%
  summarize(Mean = mean(spike),
            SD = sd(spike)) %>%
    pivot_longer(!c(delta0, sigma.delta), names_to = "stat", values_to = "count") %>%
    arrange(stat) %>%
    pivot_wider(names_from = delta0, 
                values_from = count) %>%
    relocate(stat))

print(xtable(spike.table, caption = "Spike table", type = "latex", digits = 4),
      file = "TeX/spike.table.tex", include.rownames = FALSE)
  
  
  
(spike.box <- spike %>%
    filter(delta0 != 0.5) %>%
  ggplot(aes(x = factor(delta0), y = spike)) +
  geom_boxplot() +
  facet_wrap(~ sigma.delta, nrow = 1))

