##### Summary of CTS simulation

library(knitr)
library(tidyverse)
library(xtable)
library(extrafont)
#library(gridExtra)
library(here)
library(ggplot2)
library(ggpattern)
library(scales)
library(here)
#library(simhelpers)

# cts.5 <- readRDS(here("R", "Results", "cts.sim.5.rds"))[[1]]
# cts.10 <- readRDS(here("R", "Results", "cts.sim.10.rds"))[[1]] 
# cts.30 <- readRDS(here("R", "Results", "cts.sim.30.rds"))[[1]] 
# cts.50 <- readRDS(here("R", "Results", "cts.sim.50.rds"))[[1]]
# cts.S.N.2 <- readRDS(here("R", "Results", "cts.sim.S.N.2.rds"))[[1]]
# cts.S.N.3 <- readRDS(here("R", "Results", "cts.sim.S.N.3.rds"))[[1]]

# cts.5 <- readRDS(here("R", "Results", "cts.sim.5.postprocess.rds"))[[1]] 
cts.10 <- readRDS(here("R", "Results", "cts.sim.10.postprocess.rds"))[[1]]
cts.30 <- readRDS(here("R", "Results", "cts.sim.30.postprocess.rds"))[[1]]
cts.50 <- readRDS(here("R", "Results", "cts.sim.50.postprocess.rds"))[[1]]
targets <- readRDS(here("R", "targets.rds"))
cts.all <- bind_rows(cts.10, cts.30, cts.50) %>%
  mutate_if(is.character,
            str_replace_all,
            pattern = c("sens", "spec"),
            replace = c("Sens", "Spec")) %>%
  type_convert() %>% 
  left_join(targets, by = c("CTS", "sigma.delta")) %>%
  mutate_if(is.character,
            str_replace_all,
            pattern = c("LRm", "LRp"),
            replace = c("LR-", "LR+"))


# cts.names <- c("LRm", "LRp", "NPV", "PPV", "Sens", "Spec")

# cts <- list()
# cts.S.N <- list()
# 
# for(i in 1:length(cts.10)){
# 
#   cts[[i]] <- cts.5[[i]] %>% 
#     drop_na() %>% 
#     type_convert() %>%
#     mutate(CTS = cts.names[i]) %>%
#     relocate(CTS) %>%
#     rbind(cts.10[[i]] %>% 
#             drop_na() %>% 
#             type_convert() %>%
#             mutate(CTS = cts.names[i]) %>%
#             relocate(CTS)) %>%
#     rbind(cts.30[[i]] %>% 
#             drop_na() %>% 
#             type_convert() %>%
#             mutate(CTS = cts.names[i]) %>%
#             relocate(CTS)) %>%
#     rbind(cts.50[[i]] %>% 
#             drop_na() %>% 
#             type_convert() %>%
#             mutate(CTS = cts.names[i]) %>%
#             relocate(CTS)) %>%
#     left_join(targets, by = c("CTS", "sigma.delta")) %>%
#     type_convert(guess_integer = TRUE)
#   
#   # cts.S.N[[i]] <- cts.S.N.2[[i]] %>% 
#   #   drop_na() %>% 
#   #   type_convert() %>%
#   #   mutate(CTS = cts.names[i]) %>%
#   #   relocate(CTS) %>%
#   #   rbind(cts.S.N.3[[i]] %>% 
#   #           drop_na() %>% 
#   #           type_convert() %>%
#   #           mutate(CTS = cts.names[i]) %>%
#   #           relocate(CTS))
#   
# }

rm(cts.5, cts.10, cts.30, cts.50)

# names(cts) <- names(cts.S.N) <- cts.names
# 
# all.cts <- do.call(rbind, cts)
# all.cts.S.N <- do.call(rbind, cts.S.N)

full.summary <- cts.all %>%
  group_by(CTS, S, sigma.delta, Method) %>%
  mutate(m = mean(mean.est),
         v = var(mean.est),
         est.bar_j = 1 / (n() - 1) * (n() * m - mean.est),
         S2Tj = 1 / (n() - 2) * ((n() - 1) * v - (n() / (n() - 1)) * (mean.est - m) ^ 2),
         rmse_j = sqrt((est.bar_j - target.mean)^2 + S2Tj)
  ) %>%
  summarize(Bias = mean(mean.est - target.mean),
            RMSE = sqrt(mean((mean.est - target.mean)^2)),
            V = var(mean.est),
            Avg.SD = mean(sd.est),
            `95% Coverage` = mean(ci.lower < target.mean & ci.upper > target.mean),
            `95% Length` = mean(ci.upper - ci.lower),
            mcse.Bias = sqrt(var(mean.est) / n()),
            mcse.RMSE = sqrt((n() - 1) * mean((rmse_j - RMSE) ^ 2)),
            mcse.Coverage = sqrt(`95% Coverage` * (1 - `95% Coverage`) / n()),
            mcse.Length = sqrt(var(ci.upper - ci.lower) / n()),
            lower.Bias = Bias - 1.96 * mcse.Bias,
            upper.Bias = Bias + 1.96 * mcse.Bias,
            lower.RMSE = RMSE - 1.96 * mcse.RMSE,
            upper.RMSE = RMSE + 1.96 * mcse.RMSE,
            lower.Coverage = `95% Coverage` - 1.96 * mcse.Coverage,
            upper.Coverage = `95% Coverage` + 1.96 * mcse.Coverage,
            lower.Length = `95% Length` - 1.96 * mcse.Length,
            upper.Length = `95% Length` + 1.96 * mcse.Length
    ) 

breaks_fun <- function(x){
  if(max(abs(x) > 0.1)){
    seq(-0.005, 0.005, 0.005)
  }
  else seq(-0.1, 0.1, 0.1)
}
# full.summary.S.N <- all.cts.S.N %>%
#   group_by(CTS, S, sigma.delta, Method) %>%
#   summarize(bias = mean(mean.est - target.mean),
#             rmse = sqrt(mean((mean.est - target.mean)^2)),
#             mcse = sqrt(v / 2500),
#             avg.sd = mean(sd.est),
#             cover = mean(ci.lower < target.mean & ci.upper > target.mean),
#             length = mean(ci.upper - ci.lower),
#             min = min(mean.est),
#             max = max(mean.est)) 

# yint.dummy <- data.frame(stat = c("Bias", "95% Coverage", "95% Length", "RMSE"), 
#                          Y = c(0, .95, NA, NA))
# 
# (LRp.plot <- full.summary %>%
#   filter(CTS == "LRp") %>%
#   pivot_longer(cols = c("Bias", "RMSE", "95% Coverage", "95% Length"), names_to = "stat") %>%
#   ggplot(aes(x = as.factor(S), y = value, group = interaction(Method, as.factor(sigma.delta)), 
#                color = as.factor(sigma.delta), linetype = Method)) +
#     geom_point() +
#     geom_line() + 
#     geom_hline(data = yint.dummy, aes(yintercept = Y), alpha = 0.5, linetype = 2) +
#     facet_wrap(~stat, scales = "free") +
#     theme_bw() +
#     labs(x = "No. studies",
#          y = NULL,
#          color = expression(sigma[delta])) +
#     scale_color_manual(values = c("deepskyblue2", "orchid3", "firebrick2")))
  

(bias.plot <- full.summary %>%
  ggplot(aes(x = as.factor(S), y = Bias, #group = interaction(Method, factor(sigma.delta)), 
             fill = factor(sigma.delta), 
             pattern = Method)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(), color = "black", pattern_fill = "white",
                   pattern_density = 0.4) +
  geom_errorbar(aes(ymin = lower.Bias, ymax = upper.Bias), position = position_dodge(.9), width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = 2, alpha = 0.5) +
  facet_wrap(~CTS, nrow = 2, ncol = 3,
             scales = "free") +
  labs(x = "Number of studies", y = "Bias") +
  theme_bw() +
  theme(text = element_text(size = 12, family = "LM Roman 10"),
        panel.spacing.x = unit(2, "lines")) +
  scale_pattern_manual(values = c("none", "circle"),
                       labels = c(expression(CTS[0]), expression(CTS[h]))) +
    scale_y_continuous(n.breaks = 3) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = expression(sigma[delta])) +
  coord_flip())




# (bias.plot <- full.summary %>%
#     ggplot(aes(x = as.factor(S), y = Bias, group = interaction(Method, factor(sigma.delta)), 
#                color = factor(sigma.delta), linetype = Method)) +
#     geom_point(position = position_dodge(width = 0.2)) +
#     geom_line(position = position_dodge(width = 0.2)) +
#     geom_hline(yintercept = 0) +
#     geom_linerange(aes(ymin = lower.Bias, ymax = upper.Bias), position = position_dodge(width = 0.2), linetype = 1) + 
#     facet_wrap(~CTS, nrow = 2, ncol = 3,
#                scales = "free") +
#     labs(x = "Number of studies", y = "Bias") +
#     theme_bw() +
#     theme(text = element_text(size = 12, family = "LM Roman 10")) +
#     scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = expression(sigma[delta])))

(cover.plot <- full.summary %>%
    ggplot(aes(x = factor(S), y = `95% Coverage`, group = interaction(Method, factor(sigma.delta)),
               linetype = Method, color = as.factor(sigma.delta))) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = lower.Coverage, ymax = upper.Coverage), position = position_dodge(width = 0.5), 
                  linetype = 1, width = 0.3) + 
    geom_hline(yintercept = 0.95) +
    geom_vline(xintercept = c(1.5, 2.5), linetype = 2,
               alpha = 0.25) +
    facet_wrap(~CTS, nrow = 2, ncol = 3, scales = "free") +
    ylim(c(0.7, 1)) +
    theme_bw() +
    labs(x = "Number of studies") +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    scale_linetype_manual(values = c(1, 3), 
                          labels = c(expression(CTS[0]), expression(CTS[h]))) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = expression(sigma[delta])))
  
(rmse.plot <- full.summary %>%
    ggplot(aes(x = as.factor(S), y = RMSE, group = interaction(Method, factor(sigma.delta)),
               linetype = Method, color = as.factor(sigma.delta))) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = lower.RMSE, ymax = upper.RMSE), position = position_dodge(width = 0.5), 
                  linetype = 1, width = 0.3) + 
    geom_vline(xintercept = c(1.5, 2.5), linetype = 2, alpha = 0.25) +
    facet_wrap(~CTS, nrow = 2, ncol = 3, scales = "free") +
    labs(x = "Number of studies", y = "RMSE") +
    theme_bw() +
    scale_linetype_manual(values = c(1, 3),
                          labels = c(expression(CTS[0]), expression(CTS[h]))) +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = expression(sigma[delta])) )

(length.plot <- full.summary %>%
    ggplot(aes(x = as.factor(S), y = `95% Length`, group = interaction(Method, factor(sigma.delta)),
               linetype = Method, color = as.factor(sigma.delta))) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = lower.Length, ymax = upper.Length), position = position_dodge(width = 0.5), 
                  linetype = 1, width = 0.3) + 
    geom_vline(xintercept = c(1.5, 2.5), linetype = 2,
               alpha = 0.25) +
    facet_wrap(~CTS, nrow = 2, ncol = 3, scales = "free") +
    labs(x = "Number of studies") +
    theme_bw() +
    scale_linetype_manual(values = c(1, 3),
                          labels = c(expression(CTS[0]), expression(CTS[h]))) +
    scale_y_continuous(n.breaks = 3) +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = expression(sigma[delta])))

ggsave("TeX/bias_plot.pdf", plot = bias.plot,
       width = 6.5, height = 7, units = "in")
ggsave("TeX/cover_plot.pdf", plot = cover.plot, 
       device = "pdf", width = 6, height = 6, units = "in")
ggsave("TeX/rmse_plot.pdf", plot = rmse.plot, 
       device = "pdf", width = 6, height = 6, units = "in")
ggsave("TeX/length_plot.pdf", plot = length.plot, 
       device = "pdf", width = 6, height = 6, units = "in")
