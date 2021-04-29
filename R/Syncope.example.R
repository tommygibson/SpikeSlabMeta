#### Spike/slab + Confusion numerical example with syncope data

setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta/Code")

set.seed(1234)
library(R2jags)
library(knitr)
library(gplots)
library(pander)
library(tidyverse)
library(extrafont)

univOR <- read.csv("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/Research Weiss/Read into R/For JAGS/Data/Univariate_jags2.csv")

#get rid of obs with no info
#get rid of obervation(s) from Oh's paper -- weird estimate for ECG
univOR <- filter(univOR, 
                 event_exposure > 0 | Orhat > 0,
                 !(Paper %in% c("Oh", "Colivicchi", "Del Rosso", "Martin")),
                 !(Variable %in% c("Myocardial Infarction", "Incontinence")),
                 !(Variable == "ECG" & Paper == "Quinn"),
                 !(Variable == "Hemoglobin" & Paper == "Thiruganasambandamoorthy"))

# gotta get rid of Colivicci, del rosso, martin as well

#Derose was a subset of Gabayan, and those were the only two papers to contribute to the meta analysis for myocardial infarction
# Also, only one study has incontinence


# We're taking ECG data from quinn's 2011 paper, which uses the same data as his 2004 paper
#so delete ECG info from the 2004 paper
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
                                      TRUE ~ Variable),
                 Varnum = as.numeric(factor(Variable, levels = unique(Variable))),
                 Paper = as.factor(as.character(Paper)),
                 Paper_num = as.numeric(Paper),
                 counts = is.na(event_exposure) + 1,
                 n_i0 = event_noexposure + noevent_noexposure,
                 n_i1 = event_exposure + noevent_exposure,
                 y_i0 = event_noexposure,
                 y_i1 = event_exposure,
                 N_i = n_i0 + n_i1)

keeps <- c("Variable", "lnORhat", "SE.lnOR", "n_i0", "n_i1", "y_i0", "y_i1", "N_i", "Paper_num", "Paper", "counts", "Varnum")                 
syncope <- univOR %>% 
  group_by(Variable) %>%
  arrange(Varnum, desc(lnORhat)) %>%
  select(keeps)


######################################
######### MODELS
######################################

sink("meta_confusion.txt")
cat("
model
{
  for(i in 1:S){
    y[i,1] ~ dbin(pi[i,1], n[i,1])
    y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    # meta-regression on the log-odds ratio
    # do we include intercept?
    
    delta[i] ~ dnorm(delta0, D.delta)
    
    nu[i] ~ dnorm(nu0, D.nu)
    
    beta[i] ~ dnorm(beta0, D.beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  delta0 ~ dnorm(c, d)
  nu0 ~ dnorm(e, f)
  
  # HALF T ON ALL THESE BITCHES
  D.beta <- pow(sigma.beta, -2)
  D.nu <- pow(sigma.nu, -2)
  D.delta <- pow(sigma.delta, -2)
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 15) T(0,)
  sigma.nu ~ dt(0, 1, 15) T(0,)
  sigma.delta ~ dt(0, 1, 15) T(0,)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D.delta)
    betanew[j] ~ dnorm(beta0, D.beta)
    nunew[j] ~ dnorm(nu0, D.nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi1new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))
    pi0new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))
    
    pi11new[j] <- pi1new[j] * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - pi1new[j]) * psinew[j]            # P(no event, risk)
    pi01new[j] <- pi0new[j] * (1 - psinew[j])            # P(event, no risk)
    pi00new[j] <- (1 - pi0new[j]) * (1 - psinew[j])      # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
  pi1.h <- 1 / (1 + exp(-(beta0 + delta0 / 2)))
  pi0.h <- 1 / (1 + exp(-(beta0 - delta0 / 2)))
  psi.h <- 1 / (1 + exp(-(nu0)))
  
  pi11.h <- pi1.h * psi.h
  pi10.h <- (1 - pi1.h) * psi.h
  pi01.h <- pi0.h * (1 - psi.h)
  pi00.h <- (1 - pi0.h) * (1 - psi.h)
  
  sens.h <- pi11.h / (pi11.h + pi01.h)
  spec.h <- pi00.h / (pi00.h + pi10.h)
  
  LRp.h <- sens.h / (1 - spec.h)
  LRm.h <- (1 - sens.h) / spec.h
  
}", fill = TRUE)
sink()

sink("meta_confusion_gamma.txt")
cat("
model
{
  for(i in 1:S){
    y[i,1] ~ dbin(pi[i,1], n[i,1])
    y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    # meta-regression on the log-odds ratio
    # do we include intercept?
    
    delta[i] ~ dnorm(delta0, D.delta)
    
    nu[i] ~ dnorm(nu0, D.nu)
    
    beta[i] ~ dnorm(beta0, D.beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  delta0 ~ dnorm(c, d)
  nu0 ~ dnorm(e, f)
  
  # Precisions
  D.beta ~ dgamma(3, 2)
  D.nu ~ dgamma(3, 2)
  D.delta ~ dgamma(3, 2)
  
  # gamma prior on precisions
  sigma.beta <- pow(D.beta, -1/2)
  sigma.nu <- pow(D.nu, -1/2)
  sigma.delta <- pow(D.delta, -1/2)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D.delta)
    betanew[j] ~ dnorm(beta0, D.beta)
    nunew[j] ~ dnorm(nu0, D.nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi1new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))
    pi0new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))
    
    pi11new[j] <- pi1new[j] * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - pi1new[j]) * psinew[j]            # P(no event, risk)
    pi01new[j] <- pi0new[j] * (1 - psinew[j])            # P(event, no risk)
    pi00new[j] <- (1 - pi0new[j]) * (1 - psinew[j])      # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
  pi1.h <- 1 / (1 + exp(-(beta0 + delta0 / 2)))
  pi0.h <- 1 / (1 + exp(-(beta0 - delta0 / 2)))
  psi.h <- 1 / (1 + exp(-(nu0)))
  
  pi11.h <- pi1.h * psi.h
  pi10.h <- (1 - pi1.h) * psi.h
  pi01.h <- pi0.h * (1 - psi.h)
  pi00.h <- (1 - pi0.h) * (1 - psi.h)
  
  sens.h <- pi11.h / (pi11.h + pi01.h)
  spec.h <- pi00.h / (pi00.h + pi10.h)
  
  LRp.h <- sens.h / (1 - spec.h)
  LRm.h <- (1 - sens.h) / spec.h
  
}", fill = TRUE)
sink()

sink("meta_confusion_spike.txt")
cat("
model
{
  for(i in 1:S){
    y[i,1] ~ dbin(pi[i,1], n[i,1])
    y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    
    delta[i] ~ dnorm(delta0, D.delta)
    
    nu[i] ~ dnorm(nu0, D.nu)
    
    beta[i] ~ dnorm(beta0, D.beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  nu0 ~ dnorm(e, f)
  
  delta0 <- delta1 * rho
  delta1 ~ dnorm(c, d)
  rho ~ dbern(p)
  spike <- 1 - rho
  
  # HALF T ON ALL THESE BITCHES
  D.beta <- pow(sigma.beta, -2)
  D.nu <- pow(sigma.nu, -2)
  D.delta <- pow(sigma.delta, -2)
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 1) T(0,)
  sigma.nu ~ dt(0, 1, 1) T(0,)
  sigma.delta ~ dt(0, 1, 1) T(0,)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D.delta)
    betanew[j] ~ dnorm(beta0, D.beta)
    nunew[j] ~ dnorm(nu0, D.nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi1new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))
    pi0new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))
    
    pi11new[j] <- pi1new[j] * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - pi1new[j]) * psinew[j]            # P(no event, risk)
    pi01new[j] <- pi0new[j] * (1 - psinew[j])            # P(event, no risk)
    pi00new[j] <- (1 - pi0new[j]) * (1 - psinew[j])      # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
}", fill = TRUE)
sink()

sink("meta_confusion_spike_gamma.txt")
cat("
model
{
  for(i in 1:S){
    y[i,1] ~ dbin(pi[i,1], n[i,1])
    y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    
    delta[i] ~ dnorm(delta0, D.delta)
    
    nu[i] ~ dnorm(nu0, D.nu)
    
    beta[i] ~ dnorm(beta0, D.beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  nu0 ~ dnorm(e, f)
  
  delta0 <- delta1 * rho
  delta1 ~ dnorm(c, d)
  rho ~ dbern(p)
  spike <- 1 - rho
  
  # Precisions
  D.beta ~ dgamma(3, 2)
  D.nu ~ dgamma(3, 2)
  D.delta ~ dgamma(3, 2)
  
  # gamma prior on precisions
  sigma.beta <- pow(D.beta, -1/2)
  sigma.nu <- pow(D.nu, -1/2)
  sigma.delta <- pow(D.delta, -1/2)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D.delta)
    betanew[j] ~ dnorm(beta0, D.beta)
    nunew[j] ~ dnorm(nu0, D.nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi1new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))
    pi0new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))
    
    pi11new[j] <- pi1new[j] * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - pi1new[j]) * psinew[j]            # P(no event, risk)
    pi01new[j] <- pi0new[j] * (1 - psinew[j])            # P(event, no risk)
    pi00new[j] <- (1 - pi0new[j]) * (1 - psinew[j])      # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
}", fill = TRUE)
sink()



# functions to generate initial values

init.gen <- function(){
  list(
    delta = runif(S, -1, 1),
    beta = runif(S, -1, 1),
    nu = runif(S, -1, 1),
    delta0 = runif(1, -1, 1),
    beta0 = runif(1, -1, 1),
    nu0 = runif(1, -1, 1),
    sigma.delta = runif(1, 0.5, 1),
    sigma.beta = runif(1, 0.5, 1),
    sigma.nu = runif(1, 0.5, 1),
    deltanew = runif(M, -1, 1),
    betanew = runif(M, -1, 1),
    nunew = runif(M, -1, 1)
    # delta = runif(S, -.1, .1),
    # beta = runif(S, -.1, .1),
    # nu = runif(S, -.1, .1),
    # delta0 = runif(1, -.1, .1),
    # beta0 = runif(1, -.1, .1),
    # nu0 = runif(1, -.1, .1),
    # sigma.delta = runif(1, 0.5, 1),
    # sigma.beta = runif(1, 0.5, 1),
    # sigma.nu = runif(1, 0.5, 1),
    # deltanew = rep(0, M),
    # betanew = rep(0, M),
    # nunew = rep(0, M)
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
    sigma.delta = runif(1, 0.2, 1),
    sigma.beta = runif(1, 0.2, 1),
    sigma.nu = runif(1, 0.2, 1),
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

set.seed(11223)
index <- 0
for(i in 1:max(syncope$Varnum)){

  syncope_curr <- syncope %>%
    filter(Varnum == i, counts == 1)
  
  # only do the analysis with at least two papers that provide counts
  if(dim(syncope_curr)[1] < 2) next 
  
  index <- index + 1
  
  # data for the analyses
  M <- 100
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
  
  meta.data <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                          a = -2, b = 0.5, c = 0, d = 0.5, e = 0, f = 0.5)
  
  meta.data.spike <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                    a = -2, b = 0.5, c = 0, d = 0.1, e = 0, f = 0.5, p = 0.5)

  # as a first run we'll just follow the hyperparameters
  meta.params <- c("LRmnew", "LRpnew", "PPVnew", "NPVnew", "sensnew", "specnew")
  # meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  meta.params.spike <- c("spike")
  
  #### If there are <5 studies, we use the model with informative gamma prior on precisions
  #### otherwise we use half-cauchy prior
  
  which.model <- ifelse(dim(syncope_curr)[1] < 6, 2, 1)
  # half-t prior on SDs of REs
  if(which.model == 1){
    
    meta.anal.spike <- jags(data = meta.data.spike, inits = init.gen.spike, parameters.to.save = meta.params.spike,
                            model.file = "meta_confusion_spike.txt",
                            n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
    
  }
  # gamma prior on precision of REs
  else if(which.model == 2){
    meta.anal.spike <- jags(data = meta.data.spike, inits = init.gen.spike.gamma, parameters.to.save = meta.params.spike,
                            model.file = "meta_confusion_spike_gamma.txt",
                            n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  }
  
  # look at spikes: if spike at 0 is big (> 0.5) we won't look at CTSs
  # otherwise, do full 3RE model with either gamma or half-t prior
  
  # half-t
  if(meta.anal.spike$BUGSoutput$summary[1] < 0.25 & which.model == 1){
    
    meta.anal <- jags(data = meta.data, inits = init.gen, parameters.to.save = meta.params,
                      model.file = "meta_confusion.txt",
                      n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
    meta.sims <- meta.anal$BUGSoutput$sims.matrix
    summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), make_CTS_sum(meta.sims, M)), nrow = 1)
    
    # gamma
  } else if(meta.anal.spike$BUGSoutput$summary[1] < 0.25 & which.model == 2){
    
    meta.anal <- jags(data = meta.data, inits = init.gen.gamma, parameters.to.save = meta.params,
                      model.file = "meta_confusion_gamma.txt",
                      n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
    meta.sims <- meta.anal$BUGSoutput$sims.matrix
    summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), make_CTS_sum(meta.sims, M)), nrow = 1)
    
    
  } else if(meta.anal.spike$BUGSoutput$summary[1] >= 0.25){
    
    summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], sigfig(meta.anal.spike$BUGSoutput$summary[1]), rep("--", length(meta.params))), nrow = 1)
    
  }
  
  
  colnames(summaries[[index]]) <- c("Variable", "Num.papers", "P(=0)", meta.params[order(meta.params)])
  
  # summaries.spike[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_simple_sum(meta.anal.spike$BUGSoutput$summary)), nrow = 1)
  
  
  
}

syncope_summary <- as.data.frame(do.call(rbind, summaries))[, c(1:3, 5, 4, 6:9)]
names(syncope_summary) <- c("Variable", "Num. Papers", "Spike", "LR+", "LR-", "NPV", "PPV", "Sens", "Spec")

# save results!
print(xtable(syncope_summary, caption = "Results of 31 meta-analyses of syncope studies", type = "latex"), file = "syncope.summary.tex", include.rownames = FALSE)

######## We want to test why these guys are spitting out errors

# which_errors <- c(8, 10, 16, 19, 22, 23, 25, 30, 31, 32)
# syncope_baddies <- syncope %>%
#   filter(Varnum %in% which_errors)
# 
# summaries_baddies <- list()
# 
# set.seed(11223)
# index <- 0
# for(i in 1:max(syncope$Varnum)){
#   
#   syncope_curr <- syncope %>%
#     filter(Varnum == i, counts == 1)
#   
#   # only do the analysis with at least two papers that provide counts
#   if(dim(syncope_curr)[1] < 2) next 
#   
#   index <- index + 1
#   
#   # syncope_curr <- rbind(syncope_curr, syncope_curr)
#   # data for the analyses
#   M <- 100
#   S <- dim(syncope_curr)[1]
#   n <- syncope_curr %>%
#     ungroup() %>%
#     select(c("n_i1", "n_i0"))
#   y <- syncope_curr %>%
#     ungroup() %>%
#     select(c("y_i1", "y_i0"))
#   n.tot = rowSums(n)
#   
#   
#   # we'll do two analyses per variable
#   # one with the spike/slab prior, one without
#   
#   # these values will be the same for pretty much every anal
#   # until we do sensitivity analysis
#   
#   meta.data <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
#                     a = -2, b = 0.5, c = 0, d = 0.5, e = 0, f = 0.5)
#   
#   # meta.data.spike <- list(M = M, S = S, n = n, y = y, n.tot = n.tot, 
#   #                   a = -2, b = 0.5, c = 0, d = 0.1, e = 0, f = 0.5, p = 0.5)
#   
#   # as a first run we'll just follow the hyperparameters
#   meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
#   # meta.params.spike <- c("beta0", "delta0", "nu0", "rho")
#   
#   meta.anal <- jags(data = meta.data, inits = init.gen, parameters.to.save = meta.params,
#                     model.file = "meta_confusion.txt",
#                     n.chains = 2, n.iter = 10000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
#   # meta.anal.spike <- jags(data = meta.data.spike, inits = init.gen.spike, parameters.to.save = meta.params.spike,
#   #                         model.file = "meta_confusion_spike.txt", 
#   #                         n.chains = 2, n.iter = 21000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
#   
#   # meta.sum <- meta.anal.spike$BUGSoutput$summary
#   
#   # meta.sims <- meta.anal$BUGSoutput$sims.matrix
#   
#   
#   
#   # summaries[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_CTS_sum(meta.anal$BUGSoutput$sims.matrix, 1)), nrow = 1)
#   summaries_baddies[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_simple_sum(meta.anal$BUGSoutput$summary)), nrow = 1)
#   
#   
#   colnames(summaries_baddies[[index]]) <- c("Variable", "Num.papers", meta.params[order(meta.params)])
#   
#   # summaries[[i]] <- c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_sum(meta.sum))
#   
#   
# }
# 
# do.call(rbind, summaries_baddies)


#### Try again but with gamma priors on the standard devations

summaries_gamma <- list()

set.seed(11223)
index <- 0
for(i in 1:max(syncope$Varnum)){
  
  syncope_curr <- syncope %>%
    filter(Varnum == i, counts == 1)
  
  # only do the analysis with at least two papers that provide counts
  if(dim(syncope_curr)[1] < 2) next 
  
  index <- index + 1
  
  # syncope_curr <- rbind(syncope_curr, syncope_curr)
  # data for the analyses
  M <- 100
  S <- dim(syncope_curr)[1]
  n <- syncope_curr %>%
    ungroup() %>%
    select(c("n_i1", "n_i0"))
  y <- syncope_curr %>%
    ungroup() %>%
    select(c("y_i1", "y_i0"))
  n.tot = rowSums(n)
  
  
  # we'll do two analyses per variable
  # one with the spike/slab prior, one without
  
  # these values will be the same for pretty much every anal
  # until we do sensitivity analysis
  
  meta.data <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                    a = -2, b = 0.5, c = 0, d = 0.5, e = 0, f = 0.5)
  
  # meta.data.spike <- list(M = M, S = S, n = n, y = y, n.tot = n.tot, 
  #                   a = -2, b = 0.5, c = 0, d = 0.1, e = 0, f = 0.5, p = 0.5)
  
  # as a first run we'll just follow the hyperparameters
  # meta.params <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")
  meta.params <- c("LRpnew", "LRmnew")
  # meta.params.spike <- c("beta0", "delta0", "nu0", "rho")
  
  meta.anal <- jags(data = meta.data, inits = init.gen.gamma, parameters.to.save = meta.params,
                    model.file = "meta_confusion_gamma.txt",
                    n.chains = 2, n.iter = 5000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  # meta.anal.spike <- jags(data = meta.data.spike, inits = init.gen.spike, parameters.to.save = meta.params.spike,
  #                         model.file = "meta_confusion_spike.txt", 
  #                         n.chains = 2, n.iter = 21000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
  
  # meta.sum <- meta.anal.spike$BUGSoutput$summary
  
  # meta.sims <- meta.anal$BUGSoutput$sims.matrix
  
  
  
  summaries_gamma[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_CTS_sum(meta.anal$BUGSoutput$sims.matrix, M)), nrow = 1)
  # summaries_gamma[[index]] <- matrix(c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_simple_sum(meta.anal$BUGSoutput$summary)), nrow = 1)
  
  
  colnames(summaries_gamma[[index]]) <- c("Variable", "Num.papers", meta.params[order(meta.params)])
  
  # summaries[[i]] <- c(syncope_curr$Variable[1], dim(syncope_curr)[1], make_sum(meta.sum))
  
  
}

do.call(rbind, summaries_gamma)



########### Example from Chu 2009

n1 <- c(11, 8, 4, 15, 5, 29, 16, 37, 13, 18)
n2 <- c(46, 38, 18, 45, 15, 169, 33, 235, 58, 24)

a <- c(9, 3, 3, 3, 0, 7, 12, 23, 8, 16)
d <- c(44, 32, 16, 44, 15, 167, 29, 230, 53, 22)
b <- n2 - d
c <- n1 - a

y <- cbind(a, c)
n <- cbind(a + b, c + d)

S <- length(a) 
n.tot <- rowSums(n)
M <- 100

# inits

meta.data <- list(S = S, y = y, n = n, M = M, n.tot = n.tot, a = -2, b = 0.5, c = 1, d = 0.5, 
                  e = -1, f = 0.5, p = 0.5)



meta.params <- c("")

meta_anal1 <- jags(meta.data, inits = init.gen.spike, meta.params, "meta_confusion_spike.txt", n.chains=3,
                   n.iter=5000, n.burnin=1000, n.thin=2, DIC=F)

make_CTS_sum(meta_anal1$BUGSoutput$sims.matrix[,2:201], M)
meta_anal1$BUGSoutput$summary[1:3,]

meta_sims1 <- meta_anal1$BUGSoutput$sims.matrix

cor(as.vector(sens.logits), as.vector(spec.logits))


