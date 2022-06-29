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

source(here("R", "spike.functions.R"))

univOR <- read.csv(here("R", "syncope_data.csv"))

#get rid of obs with no info
#get rid of obervation(s) from Oh's paper -- weird estimate for ECG
# gotta get rid of Colivicci, del rosso, martin as well

#Derose was a subset of Gabayan, and those were the only two papers to contribute to the meta analysis for myocardial infarction
# Also, only one study has incontinence


# We're taking ECG data from quinn's 2011 paper, which uses the same data as his 2004 paper
#so delete ECG info from the 2004 paper
univOR <- filter(univOR, 
                 event_exposure > 0 | Orhat > 0,
                 !(Paper %in% c("Oh", "Colivicchi", "Del Rosso", "Martin")),
                 !(Variable %in% c("Myocardial Infarction", "Incontinence")),
                 !(Variable == "ECG" & Paper == "Quinn"),
                 !(Variable == "Hemoglobin" & Paper == "Thiruganasambandamoorthy"))

#fill in contingency table counts
univOR <- univOR %>% mutate(
  event_noexposure = event_total - event_exposure,
  noevent_noexposure = nonevent_total - noevent_exposure,
  ORhat1 = (event_exposure * noevent_noexposure) / (event_noexposure * noevent_exposure),
  ORhat2 = ifelse(is.na(Orhat), ORhat1, Orhat),
  lnORhat = log(ORhat2),
  ln.lower = log(OR.lower),
  ln.upper = log(OR.upper),
  SE.chi = sqrt(lnORhat ^ 2 / Chisquare),
  SE.counts = sqrt(1 / event_noexposure + 1 / event_exposure + 1 / noevent_exposure + 1 / noevent_noexposure),
  SE.extrap = ((lnORhat - ln.lower) / 1.96 + (ln.upper - lnORhat) / 1.96) / 2,
  SE.lnOR = case_when(!is.na(SE.extrap) ~ SE.extrap,
                      !is.na(SE.counts) ~ SE.counts,
                      !is.na(SE.chi) ~ SE.chi),
  Variable = as.character(Variable),
  Variable = case_when(Variable == "Congestive" ~ "CHF",
                       Variable == "Hemoglobin" ~ "Hematocrit",
                       Variable == "Nonwhite Race" ~ "White Race",
                       TRUE ~ Variable),
  Varnum = as.numeric(factor(Variable, levels = unique(Variable))),
  Paper = as.factor(as.character(Paper)),
  Paper_num = as.numeric(Paper),
  counts = is.na(event_exposure) + 1,
  # reverse things for nonwhite race so that log(OR) is positive
  # i.e. it's white race as the risk factor now
  n_i0 = case_when(Variable == "White Race" ~ event_exposure + noevent_exposure,
                   Variable != "White Race" ~ event_noexposure + noevent_noexposure),
  n_i1 = case_when(Variable == "White Race" ~ event_noexposure + noevent_noexposure,
                   Variable != "White Race" ~ event_exposure + noevent_exposure),
  y_i0 = case_when(Variable == "White Race" ~ event_exposure,
                   Variable != "White Race" ~ event_noexposure),
  y_i1 = case_when(Variable == "White Race" ~ event_noexposure,
                   Variable != "White Race" ~ event_exposure),
  N_i = n_i0 + n_i1)

keeps <- c("Variable", "lnORhat", "SE.lnOR", "n_i0", "n_i1", "y_i0", "y_i1", "N_i", "Paper_num", "Paper", "counts", "Varnum")                 
syncope <- univOR %>% 
  group_by(Variable) %>%
  arrange(Varnum, desc(lnORhat)) %>%
  select(all_of(keeps))

# write.csv(syncope, 'R/syncope.cleaned.csv')



# functions to generate initial values

init.gen <- function(){
  list(
    delta = runif(S, -1, 1),
    beta = runif(S, -1, 1),
    nu = runif(S, -1, 1),
    delta0 = runif(1, -1, 1),
    beta0 = runif(1, -1, 1),
    nu0 = runif(1, -1, 1),
    sigma.delta = runif(1, 0.1, 1),
    sigma.beta = runif(1, 0.1, 1),
    sigma.nu = runif(1, 0.1, 1),
    deltanew = runif(M, -1, 1),
    betanew = runif(M, -1, 1),
    nunew = runif(M, -1, 1)
  )
}
init.gen.spike <- function(){
  list(
    delta = rnorm(S),
    beta = rnorm(S),
    nu = rnorm(S),
    delta1 = rnorm(1),
    beta0 = rnorm(1),
    nu0 = rnorm(1),
    sigma.delta = runif(1, 0.1, 1),
    sigma.beta = runif(1, 0.1, 1),
    sigma.nu = runif(1, 0.1, 1),
    deltanew = rnorm(M),
    betanew = rnorm(M),
    nunew = rnorm(M),
    rho = 1
  )
}

init.gen.gamma <- function(){
  list(
    delta = runif(S, -1, 1),
    beta = runif(S, -1, 1),
    nu = runif(S, -1, 1),
    delta0 = runif(1, -1, 1),
    beta0 = runif(1, -1, 1),
    nu0 = runif(1, -1, 1),
    D.beta = runif(1, 0.5, 1.5),
    D.delta = runif(1, 0.5, 1.5),
    D.nu = runif(1, 0.5, 1.5),
    deltanew = runif(M, -1, 1),
    betanew = runif(M, -1, 1),
    nunew = runif(M, -1, 1)
  )
}

init.gen.spike.gamma <- function(){
  list(
    delta = runif(S, -1, 1),
    beta = runif(S, -1, 1),
    nu = runif(S, -1, 1),
    delta1 = runif(1, -1, 1),
    beta0 = runif(1, -1, 1),
    nu0 = runif(1, -1, 1),
    D.beta = runif(1, 0.5, 1.5),
    D.delta = runif(1, 0.5, 1.5),
    D.nu = runif(1, 0.5, 1.5),
    deltanew = runif(M, -1, 1),
    betanew = runif(M, -1, 1),
    nunew = runif(M, -1, 1),
    rho = 1
  )
}




# trim a number to exactly n significant digits, including zeros on the end
sigfig <- function(x, n=3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
  
  trimws(format(round(x, n), nsmall = n), which = "both")
  
}   
# function to make CI from jags summary

make_ci <- function(x){
  ci <- paste("(", x[3], ", ", x[7], ")", sep = "")
  return(ci)
}


# make mean (95 %CI) for each row of output from jags summary
make_simple_sum <- function(x){
  x <- sigfig(x, 3)
  summ <- vector(length = dim(x)[1])
  for(i in 1:dim(x)[1]){
    summ[i] <- paste(c(x[i, 1], make_ci(x[i,])), collapse = " ")
  }
  return(summ)
}


make_CTS_sum <- function(x, M = 100){
  
  n.stat <- dim(x)[2] %/% M
  
  stats <- apply(x[, 1:M], 1, mean)
  if(n.stat > 1){
    for(i in 2:n.stat){
      stats <- cbind(stats, apply(x[, ((M * (i - 1)) + 1): (M * i)], 1, mean))
      
    } 
  }
  
  stats <- matrix(stats, ncol = n.stat)
  
  cond.mean <- function(x){
    return(mean(x[x < 20]))
  }
  means <- sigfig(apply(stats, 2, cond.mean), 2)
  CIs <- paste("(", apply(sigfig(apply(stats, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ", "), ")", sep = "")
  stat.summ <- apply(rbind(means, CIs), 2, paste, collapse = " ")
  
  return(stat.summ)
  
}

conditional_summary <- function(x){
  x.cond <- x[x < 30]
  return(mean(x.cond), paste)
}
make_CTS_sum <- function(x, M = 100){
  
  n.stat <- dim(x)[2] %/% M
  
  stats <- apply(x[, 1:M], 1, mean)
  if(n.stat > 1){
    for(i in 2:n.stat){
      stats <- cbind(stats, apply(x[, ((M * (i - 1)) + 1): (M * i)], 1, mean))
      
    } 
  }
  
  stats <- matrix(stats, ncol = n.stat)
  
  cond.mean <- function(x){
    return(mean(x[x < 20]))
  }
  means <- sigfig(apply(stats, 2, cond.mean), 2)
  CIs <- paste("(", apply(sigfig(apply(stats, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ", "), ")", sep = "")
  stat.summ <- apply(rbind(means, CIs), 2, paste, collapse = " ")
  
  return(stat.summ)
  
}

make_CTS_sum1 <- function(x, M = 100){
  
  # averaging over stat_new in each iteration to get the posterior
  stats <- apply(x[, 1:M], 1, mean)
  
  
  means <- mean(stats)
  low.up <- quantile(stats, c(.025, .5, .975))
  sds <- sd(stats)
  stat.summ <- c(means, sds, low.up)
  
  return(stat.summ)
  
}


# initialize matrices to hold values we care about

summaries <- list()
hyper.summaries <- list()

set.seed(11223)
index <- 0

# sims.for.densities <- list()
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
  
  # First we do spike/slab model, and if posterior spike is small then we do regular 3RE model
  
  # these values will be the same for pretty much every anal
  # until we do sensitivity analysis
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
  
  # as a first run we'll just follow the hyperparameters
  #meta.params <- c("LRmnew", "LRpnew", "PPVnew", "NPVnew", "sensnew", "specnew")
  meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  meta.params.spike <- c("spike")
  
  #### If there are <5 studies, we use the model with informative gamma prior on precisions
  #### otherwise we use half-cauchy prior
  
  # which.model <- ifelse(dim(syncope_curr)[1] < 6, 2, 1)
  # half-t prior on SDs of REs
  #if(which.model == 1){
  
  # fit 3RE-SAS and 3RE models
  meta.anal.spike <- do.call(jags.parallel,
                             list(data = names(meta.data.spike), inits = init.gen.spike, parameters.to.save = meta.params.spike,
                                  model.file = here("R", "Models", "meta_confusion_spike.txt"),
                                  n.chains = 4, n.iter = 5000, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  meta.anal <- do.call(jags.parallel,
                       list(data = names(meta.data), inits = init.gen, parameters.to.save = meta.params,
                            model.file = here("R", "Models", "meta_confusion.txt"),
                            n.chains = 4, n.iter = 5000, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  
  # save summaries of hyperparameters
  hyper.summaries[[index]] <- matrix(apply(meta.anal$BUGSoutput$summary, 1, function(x){
    paste(sigfig(x[1]), " (", sigfig(x[3]), ", ", sigfig(x[7]), ")", sep = "")
  }), nrow = 1)
  colnames(hyper.summaries[[index]]) <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  
  
  # extract posterior draws, post-processing for unreasonable values of heterogeneity parameters
  meta.sims <- as.data.frame(meta.anal$BUGSoutput$sims.matrix) %>%
    filter(sigma.beta < 5, sigma.delta < 5, sigma.nu < 5)
  
  # calculate CTSs for each combination of hyperparameters
  cts0.full <- make_cts0_from_hyper(meta.sims, M = M)
  
  # save LR+ and spec draws
  # sims.for.densities[[i]] <- cts0.full[,c(1, 2)]
  
  # posterior summaries for all CTSs
  cts0.summ <- apply(cts0.full, 2, function(x){
    paste(sigfig(mean(x[x < 30]), 2), " (", sigfig(quantile(x[x < 30], .025), 2), ", ", sigfig(quantile(x[x < 30], .975), 2), ")", sep = "")
  })
  
  # save all relevant information (including variable name, number of papers, spike height, and mean (95% CI) for CTSs)
  summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], meta.anal.spike$BUGSoutput$summary[1], cts0.summ), nrow = 1)#make_CTS_sum(meta.sims, M)), nrow = 1)
  # foo <- make_cts0_from_hyper(meta.sims, M = M)
  #}
  # gamma prior on precision of REs
  # else if(which.model == 2){
  #   meta.anal.spike <- do.call(jags.parallel,
  #                              list(data = names(meta.data.spike), inits = init.gen.spike.gamma, parameters.to.save = meta.params.spike,
  #                              model.file = here("R", "meta_confusion_spike_gamma.txt"),
  #                              n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  #   meta.anal <- do.call(jags.parallel,
  #                        list(data = names(meta.data), inits = init.gen.gamma, parameters.to.save = meta.params,
  #                        model.file = here("R", "meta_confusion_gamma.txt"),
  #                        n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  #   meta.sims <- meta.anal$BUGSoutput$sims.matrix
  #   summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), make_CTS_sum(meta.sims, M)), nrow = 1)
  # }
  
  
  # half-t
  # this commented bit is deprecated from when we only continued with analysis if the spike was small
  # if(meta.anal.spike$BUGSoutput$summary[1] < 0.25 & which.model == 1){
  #   meta.anal <- jags(data = meta.data, inits = init.gen, parameters.to.save = meta.params,
  #                     model.file = "meta_confusion.txt",
  #                     n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  #   meta.sims <- meta.anal$BUGSoutput$sims.matrix
  #   summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), make_CTS_sum(meta.sims, M)), nrow = 1)
  #   
  #   # gamma
  # } #else if(meta.anal.spike$BUGSoutput$summary[1] < 0.25 & which.model == 2){
  #   meta.anal <- jags(data = meta.data, inits = init.gen.gamma, parameters.to.save = meta.params,
  #                     model.file = "meta_confusion_gamma.txt",
  #                     n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  #   meta.sims <- meta.anal$BUGSoutput$sims.matrix
  #   summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), make_CTS_sum(meta.sims, M)), nrow = 1)
  #   
  #   
  # }
  # #  else if(meta.anal.spike$BUGSoutput$summary[1] >= 0.25){
  #   
  #   summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), rep("--", length(meta.params))), nrow = 1)
  #   
  # }
  
  
  colnames(summaries[[index]]) <- c("Variable", "Num.papers", "P(=0)", "LR+", "LR-", "NPV", "PPV", "Sens", "Spec")#meta.params[order(meta.params)])
  
  # summaries.spike[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_simple_sum(meta.anal.spike$BUGSoutput$summary)), nrow = 1)
  
  
  
}

syncope_summary <- as.data.frame(do.call(rbind, summaries))[, c(1:3, 5, 4, 6:9)]
names(syncope_summary) <- c("Variable", "Num. Papers", "Spike", "LR+", "LR-", "NPV", "PPV", "Sens", "Spec")

syncope_summary_typecorrect <- syncope_summary %>%
  type_convert(guess_integer = TRUE) %>%
  arrange(Spike)

# save results!
print(xtable(syncope_summary_typecorrect, caption = "Results of 31 meta-analyses of syncope studies", type = "latex", digits = 3), 
      file = "TeX/syncope.summary.tex", include.rownames = FALSE)


### density plots for top 10 risk factors (by spike height)

# top4 <- order(as.numeric(syncope_summary$Spike))[1:4]


# top4.sims <- sims.for.densities[top4]
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
  
  # First we do spike/slab model, and if posterior spike is small then we do regular 3RE model
  
  # these values will be the same for pretty much every anal
  # until we do sensitivity analysis
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
                            n.chains = 4, n.iter = 10000, n.thin = 2, n.burnin = 1000, DIC = FALSE))
  
  # extract posterior draws, post-processing for unreasonable values of heterogeneity parameters
  meta.sims <- as.data.frame(meta.anal$BUGSoutput$sims.matrix) %>%
    filter(sigma.beta < 5, sigma.delta < 5, sigma.nu < 5)
  
  # calculate CTSs for each combination of hyperparameters
  cts0.full <- make_cts0_from_hyper(meta.sims, M = M)
  
  # save LR+ and spec draws
  top4.sims[[index]] <- cbind.data.frame(cts0.full[,c(1, 2)], 
                                         syncope_curr$Variable[1],
                                         dim(syncope_curr)[1])
  
  index <- index + 1
  
}

density.plot.dat <- do.call(rbind, top4.sims) %>%
  rename(`LR+[0]` = LRp.new,
         `LR-[0]` = LRm.new,
         `Risk Factor` = `syncope_curr$Variable[1]`,
         `No. Studies` = `dim(syncope_curr)[1]`) %>%
  filter(`LR+[0]` < 30) #%>%


(age.contour <- density.plot.dat %>%
    filter(`Risk Factor` == "Age") %>% 
    ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
    stat_density2d(h = 0.11, n = 500, bins = 10) +
    theme_bw() +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    labs(title = "Old Age",
         x = NULL, y = NULL))


(male.contour <- density.plot.dat %>%
    filter(`Risk Factor` == "Male Gender") %>% 
    ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
    stat_density2d(h = 0.04, n = 500, bins = 10) +
    theme_bw() +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    labs(title = "Male Gender",
         x = NULL, y = NULL))

(chf.contour <- density.plot.dat %>%
    filter(`Risk Factor` == "CHF") %>% 
    ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
    geom_density_2d(h = 0.16, n = 500, bins = 10) +
    theme_bw() +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    labs(title = "CHF",
         x = NULL, y = NULL))

(heart.contour <- density.plot.dat %>%
    filter(`Risk Factor` == "Heart Disease") %>% 
    ggplot(aes(x = `LR+[0]`, y = `LR-[0]`)) + 
    geom_density_2d(h = 0.13, n = 500, bins = 10) +
    theme_bw() +
    theme(text = element_text(size = 12, family = "LM Roman 10")) +
    labs(title = "Heart disease",
         x = NULL, y = NULL))

xlab.contour <- textGrob(expression("LR+"[0]), gp = gpar(fontfamily = "LM Roman 10", fontsize = 16))
ylab.contour <- textGrob(expression("LR-"[0]), gp = gpar(fontfamily = "LM Roman 10", fontsize = 16), rot = 90)
contour.plot <- grid.arrange(age.contour, male.contour, chf.contour, heart.contour, 
                             bottom = xlab.contour, left = ylab.contour, nrow = 2)

ggsave(here("TeX", "syncope_contour.pdf"), contour.plot, 
       height = 6, width = 6, units = "in", dpi = 600)


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

# First we do spike/slab model, and if posterior spike is small then we do regular 3RE model

# these values will be the same for pretty much every anal
# until we do sensitivity analysis
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

print(xtable(trop.rf, caption = "Posterior summaries for CTSs given known information about P(RF) for a new study", type = "latex", digits = 2),
      file = "TeX/troponin.rf.tex", include.rownames = FALSE)



