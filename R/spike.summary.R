#### taking small results from cts simulation to see how many iterations are needed

library(tidyverse)
library(here)

spike <- readRDS(here("R", "Results", "meta.spike.results.rds"))[[1]]

(spike.table <- spike %>% 
  filter(delta0 != 0.5) %>%
  group_by(delta0, sigma.delta) %>%
  summarize(Mean = mean(spike),
            SD = sd(spike)))
  
  
  
(spike.box <- spike %>%
    filter(delta0 != 0.5) %>%
  ggplot(aes(x = factor(delta0), y = spike)) +
  geom_boxplot() +
  facet_wrap(~ sigma.delta, nrow = 1))

