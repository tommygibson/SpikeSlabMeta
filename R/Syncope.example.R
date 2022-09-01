#### Spike/slab + Confusion numerical example with syncope data

# setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta/Code")


set.seed(1234)
library(here)
library(R2jags)
library(knitr)
library(tidyverse)
library(extrafont)
library(xtable)
library(gridExtra)
library(grid)
library(ks)


source(here("R", "spike.functions.R"))

syncope <- read.csv(here("R", "syncope.cleaned.csv"))


# initialize matrices to hold values we care about

summaries <- list()
hyper.summaries <- list()

set.seed(11223)
index <- 0
num.removed <- list()
frac.removed <- list()
final.samples <- list()

for(i in 1:max(syncope$Varnum)){
  
  syncope_curr <- syncope %>%
    filter(Varnum == i, counts == 1)
  
  # only do the analysis with at least two papers that provide counts
  if(dim(syncope_curr)[1] < 2) next 
  
  index <- index + 1
  
  # data for the analyses
  M <- 1000
  S <- dim(syncope_curr)[1]
  n <- syncope_curr %>%
    ungroup() %>%
    select(c("n_i1", "n_i0"))
  y <- syncope_curr %>%
    ungroup() %>%
    select(c("y_i1", "y_i0"))
  n.tot = rowSums(n)
  
  # these values will be the same for pretty much every analysis
  # until we do sensitivity analysis
  a <- -2
  b <- 0.25
  c <- 0
  d <- 0.25
  e <- 0
  f <- 0.25
  p <- 0.5
  
  #### If there are <5 studies, precision on half-cauchy for heterogeneity parameters is 16 (scale = 0.25)
  #### otherwise precision is 4 (scale = 0.5)
  if(S < 5){
    B.beta <- B.nu <- B.delta <- 16
  } else {
    B.beta <- B.nu <- B.delta <- 4
  }
  
  meta.data <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                    a = a, b = b, c = c, d = d, e = e, f = f,
                    B.beta = B.beta, B.nu = B.nu, B.delta = B.delta)
  
  meta.data.spike <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                          a = a, b = b, c = c, d = d, e = e, f = f, p = p,
                          B.beta = B.beta, B.nu = B.nu, B.delta = B.delta)
  
  # follow hyperparameters in 3RE model
  # follow spike in 3RE-SAS model
  meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  meta.params.spike <- c("spike")
  
  
  
  # fit 3RE-SAS and 3RE models
  meta.anal.spike <- do.call(jags.parallel,
                             list(data = names(meta.data.spike), inits = init.gen.spike, parameters.to.save = meta.params.spike,
                                  model.file = here("R", "Models", "meta_confusion_spike.txt"),
                                  n.chains = 4, n.iter = 5000, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  meta.anal <- do.call(jags.parallel,
                       list(data = names(meta.data), inits = init.gen, parameters.to.save = meta.params,
                            model.file = here("R", "Models", "meta_confusion.txt"),
                            n.chains = 4, n.iter = 5100, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  
  # save summaries of hyperparameters
  hyper.summaries[[index]] <- matrix(apply(meta.anal$BUGSoutput$summary, 1, function(x){
    paste(sigfig(x[1]), " (", sigfig(x[3]), ", ", sigfig(x[7]), ")", sep = "")
  }), nrow = 1)
  colnames(hyper.summaries[[index]]) <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  
  
  # extract posterior draws, post-processing for unreasonable values of heterogeneity parameters
  meta.sims <- as.data.frame(meta.anal$BUGSoutput$sims.matrix)
  
  # calculate CTSs for each combination of hyperparameters
  cts0.full <- make_cts0_from_hyper(meta.sims, M = M)
  cts0.trim <- cts0.full[cts0.full[,2] < 30,]
  
  num.removed[[index]] <- nrow(cts0.full) - nrow(cts0.trim)
  frac.removed[[index]] <- num.removed[[index]] / nrow(cts0.full)
  final.samples[[index]] <- nrow(cts0.trim)
  
  
  # posterior summaries for all CTSs
  
  cts0.summ <- apply(cts0.trim, 2, function(x){
    paste(sigfig(mean(x), 2), " (", sigfig(quantile(x, .025), 2), ", ", sigfig(quantile(x, .975), 2), ")", sep = "")
  })
  
  # save all relevant information (including variable name, number of papers, spike height, and mean (95% CI) for CTSs)
  summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], meta.anal.spike$BUGSoutput$summary[1], cts0.summ), nrow = 1)
  
  
  colnames(summaries[[index]]) <- c("Variable", "Num.papers", "P(=0)", "LR+", "LR-", "NPV", "PPV", "Sens", "Spec")
  
  
  
  
}

syncope_summary <- as.data.frame(do.call(rbind, summaries))[, c(1:3, 5, 4, 6:9)]
names(syncope_summary) <- c("Variable", "Num. Papers", "Spike", "LR+", "LR-", "NPV", "PPV", "Sens", "Spec")

syncope_summary_typecorrect <- syncope_summary %>%
  type_convert(guess_integer = TRUE) %>%
  arrange(Spike)

# save results
print(xtable(syncope_summary_typecorrect, caption = "Results of 31 meta-analyses of syncope studies", type = "latex", digits = 3), 
      file = "TeX/syncope.summary.tex", include.rownames = FALSE)


### density plots for top 4 risk factors (by spike height)

# top 4 vars are Age, Male gender, CHF, Heart disease
top4 <- unique(syncope$Varnum[syncope$Variable %in% c("Age", "Male Gender", "CHF", "Heart Disease")])
top4.sims <- list()
index <- 1
for(i in top4){
  
  syncope_curr <- syncope %>%
    filter(Varnum == i, counts == 1)
  
  # only do the analysis with at least two papers that provide counts
  if(dim(syncope_curr)[1] < 2) next 
  
  # data for the analyses
  M <- 1000
  S <- dim(syncope_curr)[1]
  n <- syncope_curr %>%
    ungroup() %>%
    select(c("n_i1", "n_i0"))
  y <- syncope_curr %>%
    ungroup() %>%
    select(c("y_i1", "y_i0"))
  n.tot = rowSums(n)
  
  a <- -2
  b <- 0.25
  c <- 0
  d <- 0.25
  e <- 0
  f <- 0.25
  p <- 0.5
  
  if(S < 5){
    B.beta <- B.nu <- B.delta <- 16
  } else {
    B.beta <- B.nu <- B.delta <- 4
  }
  
  meta.data <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                    a = a, b = b, c = c, d = d, e = e, f = f,
                    B.beta = B.beta, B.nu = B.nu, B.delta = B.delta)
  
  meta.data.spike <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                          a = a, b = b, c = c, d = d, e = e, f = f, p = p,
                          B.beta = B.beta, B.nu = B.nu, B.delta = B.delta)
  
  meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  meta.params.spike <- c("spike")
  
  
  
  # and 3RE models with lots more samples for contour plots
  
  meta.anal <- do.call(jags.parallel,
                       list(data = names(meta.data), inits = init.gen, parameters.to.save = meta.params,
                            model.file = here("R", "Models", "meta_confusion.txt"),
                            n.chains = 4, n.iter = 105000, n.thin = 2, n.burnin = 5000, DIC = FALSE))
  
  # extract posterior draws, post-processing for unreasonable values of heterogeneity parameters
  meta.sims <- as.data.frame(meta.anal$BUGSoutput$sims.matrix)
  
  # calculate CTSs for each combination of hyperparameters
  cts0.full <- make_cts0_from_hyper(meta.sims, M = M)
  
  # save LR+ and LR- draws
  top4.sims[[index]] <- cbind.data.frame(cts0.full[,c(1, 2)], 
                                         syncope_curr$Variable[1],
                                         dim(syncope_curr)[1])
  
  index <- index + 1
  
}

# consolidate posterior draws for the four risk factors, rename some variables
density.plot.dat <- do.call(rbind, top4.sims) %>%
  rename(`LR+[0]` = LRp.new,
         `LR-[0]` = LRm.new,
         `Risk Factor` = `syncope_curr$Variable[1]`,
         `No. Studies` = `dim(syncope_curr)[1]`) %>%
  filter(`LR+[0]` < 30) 


#### ks package for kernel smoothing and contour plots

age.H <- 4 * Hpi(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`))
male.H <- 4 * Hpi(x = density.plot.dat %>% filter(`Risk Factor` == "Male Gender") %>% select(`LR+[0]`, `LR-[0]`))
chf.H <- 4 * Hpi(x = density.plot.dat %>% filter(`Risk Factor` == "CHF") %>% select(`LR+[0]`, `LR-[0]`))
heart.H <- 4 * Hpi(x = density.plot.dat %>% filter(`Risk Factor` == "Heart Disease") %>% select(`LR+[0]`, `LR-[0]`))
# age.H.diag <- Hpi.diag(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`))
# age.H.arb <- matrix(c(.5, 0, 0, .5), nrow = 2)
# age.H.scv <- Hscv(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`))
age.ks <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`),
              H = age.H, gridsize = c(500, 500))
male.ks <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "Male Gender") %>% select(`LR+[0]`, `LR-[0]`),
              H = male.H, gridsize = c(500, 500))
chf.ks <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "CHF") %>% select(`LR+[0]`, `LR-[0]`),
              H = chf.H, gridsize = c(500, 500))
heart.ks <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "Heart Disease") %>% select(`LR+[0]`, `LR-[0]`),
              H = heart.H, gridsize = c(500, 500))# age.ks.diag <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`),
#                    H = age.H.diag, gridsize = c(500, 500))
# age.ks.arb <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`),
#                   H = age.H.arb, gridsize = c(500, 500))
# age.ks.scv <- kde(x = density.plot.dat %>% filter(`Risk Factor` == "Age") %>% select(`LR+[0]`, `LR-[0]`),
#                   H = age.H.scv, gridsize = c(500, 500))

# pdf(file = here("TeX", "contour1.pdf"),
#     width = 7,
#     height = 7)

setEPS()
cairo_ps(file = "contour.eps", height = 7, width = 7)
par(mfrow = c(2, 2), 
    mar = c(1, 2, 3, 2),
    oma = c(5, 4, 0, 0),
    font.main = 1)
plot(age.ks, cont = c(5, 25, 50, 75, 95), col = "black",
     xlim = c(1.45, 2.9), ylim = c(0.24, 0.7), lwd = 1.5,
     xlab = NULL, ylab = NULL, family = "LM Roman 10",
     drawlabels = FALSE)
title(main = "Old Age", family = "LM Roman 10", adj = 0, line = 0.5,
      cex.main = 1.4)

plot(male.ks, cont = c(5, 25, 50, 75, 95), col = "black",
     xlim = c(1.25, 1.55), ylim = c(.62, .84), lwd = 1.5,
     xlab = NULL, ylab = NULL,
     drawlabels = FALSE, family = "LM Roman 10", xaxt = "n")
title(main = "Male Gender", family = "LM Roman 10", adj = 0, line = 0.5,
      cex.main = 1.4)
axis(1, at = c(1.3, 1.4, 1.5), family = "LM Roman 10")

plot(chf.ks, cont = c(5, 25, 50, 75, 95), col = "black",
     xlim = c(2, 5), ylim = c(0.7, .95), lwd = 1.5,
     xlab = NULL, ylab = NULL, family = "LM Roman 10", yaxt = "n",
     drawlabels = FALSE)
axis(2, at = c(.75, .8, .85, .9, .95), family = "LM Roman 10")
title(main = "CHF", family = "LM Roman 10", adj = 0, line = 0.5,
      cex.main = 1.4)

plot(heart.ks, cont = c(5, 25, 50, 75, 95), col = "black",
     xlim = c(1.45, 3.1), ylim = c(0.65, .94), lwd = 1.5,
     xlab = NULL, ylab = NULL, family = "LM Roman 10", yaxt = "n")#,
     # drawlabels = FALSE)
title(main = "Heart Failure", family = "LM Roman 10", adj = 0, line = 0.5,
      cex.main = 1.4)
axis(2, at = c(.65, .75, .85, .95), family = "LM Roman 10")
title(xlab = expression("LR+"[0]), outer = TRUE, cex.lab = 1.8, line = 2,
      family = "LM Roman 10", font.lab = 2)
title(ylab = expression("LR-"[0]), outer = TRUE, cex.lab = 1.8, line = 0.5,
      family = "LM Roman 10")
# title(xlab = expression("LR+"[0]),
#       ylab = expression("LR-"[0]),
#       cex.lab = 1.5,
#       outer = TRUE,
#       family = "LM Roman 10")
dev.off()


#### The following contour code is deprecated
# title(main = "Old Age", family = "LM Roman 10", line = 0.4, adj = 0)
# # contour plots for each risk factor
# 
# (age.contour <- density.plot.dat %>%
#     filter(`Risk Factor` == "Age") %>% 
#     ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
#     stat_density2d(h = 0.11, n = 500, bins = 10) +
#     theme_bw() +
#     theme(text = element_text(size = 12, family = "LM Roman 10")) +
#     labs(title = "Old Age",
#          x = NULL, y = NULL))
# 
# 
# (male.contour <- density.plot.dat %>%
#     filter(`Risk Factor` == "Male Gender") %>% 
#     ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
#     stat_density2d(h = 0.04, n = 500, bins = 10) +
#     theme_bw() +
#     theme(text = element_text(size = 12, family = "LM Roman 10")) +
#     labs(title = "Male Gender",
#          x = NULL, y = NULL))
# 
# (chf.contour <- density.plot.dat %>%
#     filter(`Risk Factor` == "CHF") %>% 
#     ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
#     geom_density_2d(h = 0.16, n = 500, bins = 10) +
#     theme_bw() +
#     theme(text = element_text(size = 12, family = "LM Roman 10")) +
#     labs(title = "CHF",
#          x = NULL, y = NULL))
# 
# (heart.contour <- density.plot.dat %>%
#     filter(`Risk Factor` == "Heart Disease") %>% 
#     ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
#     geom_density_2d(h = 0.13, n = 500, bins = 10) +
#     theme_bw() +
#     theme(text = element_text(size = 12, family = "LM Roman 10")) +
#     labs(title = "Heart disease",
#          x = NULL, y = NULL))
# 
# # x and y axis labels for grid.arrange function
# xlab.contour <- textGrob(expression("LR+"[0]), gp = gpar(fontfamily = "LM Roman 10", fontsize = 16))
# ylab.contour <- textGrob(expression("LR-"[0]), gp = gpar(fontfamily = "LM Roman 10", fontsize = 16), rot = 90)
# contour.plot <- grid.arrange(age.contour, male.contour, chf.contour, heart.contour, 
#                              bottom = xlab.contour, left = ylab.contour, nrow = 2)
# 
# ggsave(here("TeX", "syncope_contour.pdf"), contour.plot, 
#        height = 6, width = 6, units = "in", dpi = 600)


#### Example using Troponin:
#### we know P(RF) = 0.05, 0.1, or 0.25

troponin <- syncope %>% filter(Variable == "Troponin_99")

# only do the analysis with at least two papers that provide counts

# data for the analyses
M <- 1000
S <- dim(troponin)[1]
n <- troponin %>%
  ungroup() %>%
  select(c("n_i1", "n_i0"))
y <- troponin %>%
  ungroup() %>%
  select(c("y_i1", "y_i0"))
n.tot = rowSums(n)

a <- -2
b <- 0.25
c <- 0
d <- 0.25
e <- 0
f <- 0.25
p <- 0.5


B.beta <- B.nu <- B.delta <- 16

meta.data <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                  a = a, b = b, c = c, d = d, e = e, f = f,
                  B.beta = B.beta, B.nu = B.nu, B.delta = B.delta)


meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")


meta.anal <- do.call(jags.parallel,
                     list(data = names(meta.data), inits = init.gen, parameters.to.save = meta.params,
                          model.file = here("R", "Models", "meta_confusion.txt"),
                          n.chains = 4, n.iter = 20000, n.thin = 2, n.burnin = 10000, DIC = FALSE))

troponin.sims <- as.data.frame(meta.anal$BUGSoutput$sims.matrix) %>%
  filter(sigma.beta < 5, sigma.delta < 5)

# functions to calculate predictive CTS0 with known values of P(RF)

cts.known.rf <- function(gamma, rf.prob, M){
  beta0 <- gamma[1]
  delta0 <- gamma[2]
  
  sigma.beta <- gamma[4]
  sigma.delta <- gamma[5]
  
  beta.new <- rnorm(M, mean = beta0, sd = sigma.beta)
  delta.new <- rnorm(M, mean = delta0, sd = sigma.delta)
  
  pi1.new <- 1 / (1 + exp(-(beta.new + delta.new / 2)))
  pi0.new <- 1 / (1 + exp(-(beta.new - delta.new / 2)))
  psi <- rf.prob
  
  sens.new <- (pi1.new * psi) / (pi1.new * psi + pi0.new * ( 1 - psi))
  spec.new <- ((1 - pi0.new) * (1 - psi)) / ((1 - pi0.new) * (1 - psi) + (1 - pi1.new) * psi)
  LRm.new <- (1 - sens.new) / spec.new
  LRp.new <- sens.new / (1 - spec.new)
  
  mc.est <- apply(cbind(LRm.new, LRp.new, 1 - pi0.new, pi1.new, sens.new, spec.new), MARGIN = 2, mean)
  return(mc.est)
  
}

cts.known.rf.full <- function(gamma.matrix, rf.prob, M){
  t(apply(gamma.matrix, MARGIN = 1, cts.known.rf, rf.prob = rf.prob, M = M))
}

# calculate CTSs for known values of P(RF) \in {0.05, 0.10, 0.25}
trop.known.rf.05 <- cts.known.rf.full(troponin.sims, rf.prob = 0.05, M = 1000)
trop.known.rf.05.summary <- apply(trop.known.rf.05, 2, function(x){
  paste(sigfig(mean(x[x < 30]), 2), " (", sigfig(quantile(x[x < 30], .025), 2), ", ", sigfig(quantile(x[x < 30], .975), 2), ")", sep = "")
})
trop.known.rf.1 <- cts.known.rf.full(troponin.sims, rf.prob = 0.1, M = 1000)
trop.known.rf.1.summary <- apply(trop.known.rf.1, 2, function(x){
  paste(sigfig(mean(x[x < 30]), 2), " (", sigfig(quantile(x[x < 30], .025), 2), ", ", sigfig(quantile(x[x < 30], .975), 2), ")", sep = "")
})
trop.known.rf.25 <- cts.known.rf.full(troponin.sims, rf.prob = 0.25, M = 1000)
trop.known.rf.25.summary <- apply(trop.known.rf.25, 2, function(x){
  paste(sigfig(mean(x[x < 30]), 2), " (", sigfig(quantile(x[x < 30], .025), 2), ", ", sigfig(quantile(x[x < 30], .975), 2), ")", sep = "")
})

trop.rf <- cbind.data.frame(c(.05, .10, .25), rbind(trop.known.rf.05.summary,
                                                    trop.known.rf.1.summary,
                                                    trop.known.rf.25.summary))
rownames(trop.rf) <- NULL
names(trop.rf) <- c("P(RF)", "LR-", "LR+", "NPV", "PPV", "Sens", "Spec")

# save table of results
print(xtable(trop.rf, caption = "Posterior summaries for CTSs given known information about P(RF) for a new study", type = "latex", digits = 2),
      file = "TeX/troponin.rf.tex", include.rownames = FALSE)



