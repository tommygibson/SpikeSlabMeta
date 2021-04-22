#### Comparison of our trivariate model with chu's on real data

setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/My Papers/Spike Slab Meta/Code")

library(R2jags)
library(knitr)
library(ggplot2)

# data from Schiedler

n1 <- c(11, 8, 4, 15, 5, 29, 16, 37, 13, 18)
n2 <- c(46, 38, 18, 45, 15, 169, 33, 235, 58, 24)

a <- c(9, 3, 3, 3, 0, 7, 12, 23, 8, 16)
d <- c(44, 32, 16, 44, 15, 167, 29, 230, 53, 22)
b <- n2 - d
c <- n1 - a

Y <- cbind(a, c)
n <- cbind(a + b, c + d)

### Meta-analysis model

sink("meta_confusion.txt")
cat("
model
{
  for(i in 1:S){
    Y[i,1] ~ dbin(pi[i,1], n[i,1])
    Y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    # meta-regression on the log-odds ratio
    # do we include intercept?
    
    delta[i] ~ dnorm(delta0, D_delta)
    
    nu[i] ~ dnorm(nu0, D_nu)
    
    beta[i] ~ dnorm(beta0, D_beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  delta0 ~ dnorm(c, d)
  nu0 ~ dnorm(e, f)
  
  # HALF T ON ALL THESE BITCHES
  D_beta <- pow(sigma.beta, -2)
  D_nu <- pow(sigma.nu, -2)
  D_delta <- pow(sigma.delta, -2)
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 1) T(0,)
  sigma.nu ~ dt(0, 1, 1) T(0,)
  sigma.delta ~ dt(0, 1, 1) T(0,)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D_delta)
    betanew[j] ~ dnorm(beta0, D_beta)
    nunew[j] ~ dnorm(nu0, D_nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi11new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j]/2))) * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))) * psinew[j]          # P(no event, risk)
    pi01new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2))) * (1 - psinew[j])          # P(event, no risk)
    pi00new[j] <- (1 - 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))) * (1 - psinew[j])    # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
}", fill = TRUE)
sink()

sink("meta_confusion_diag.txt")
cat("
model
{
  for(i in 1:S){
    Y[i,1] ~ dbin(pi[i,1], n[i,1])
    Y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    # meta-regression on the log-odds ratio
    # do we include intercept?
    
    delta[i] ~ dnorm(delta0, D_delta)
    
    nu[i] ~ dnorm(nu0, D_nu)
    
    beta[i] ~ dnorm(beta0, D_beta)
    
  }
  
  beta0 ~ dnorm(a, b)
  delta0 ~ dnorm(c, d)
  nu0 ~ dnorm(e, f)
  
  # HALF T ON ALL THESE BITCHES
  D_beta <- pow(sigma.beta, -2)
  D_nu <- pow(sigma.nu, -2)
  D_delta <- pow(sigma.delta, -2)
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 1) T(0,)
  sigma.nu ~ dt(0, 1, 1) T(0,)
  sigma.delta ~ dt(0, 1, 1) T(0,)

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    deltanew[j] ~ dnorm(delta0, D_delta)
    betanew[j] ~ dnorm(beta0, D_beta)
    nunew[j] ~ dnorm(nu0, D_nu)
    psinew[j] <- 1 / (1  + exp(-nunew[j]))
    
    pi11new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j]/2))) * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))) * psinew[j]          # P(no event, risk)
    pi01new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2))) * (1 - psinew[j])          # P(event, no risk)
    pi00new[j] <- (1 - 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))) * (1 - psinew[j])    # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
}", fill = TRUE)
sink()

sink("meta_confusion_multivariate.txt")
cat("
model
{
  for(i in 1:S){
    Y[i,1] ~ dbin(pi[i,1], n[i,1])
    Y[i,2] ~ dbin(pi[i,2], n[i,2])
    
    # probability of exposure
    n[i,1] ~ dbin(psi[i], n.tot[i])
    
    psi[i] <- 1 / (1 + exp(-nu[i]))
    
    # conditional prob of event given exposure
    logit(pi[i,1]) <- beta[i] + delta[i]/2
    logit(pi[i,2]) <- beta[i] - delta[i]/2
    
    # meta-regression on the log-odds ratio
    # do we include intercept?
    
    beta[i] <- theta[i, 1]
    delta[i] <- theta[i, 2]
    nu[i] <- theta[i, 3]
    
    theta[i, 1:3] ~ dmnorm.vcov(mu[1:3], Sigma[1:3, 1:3])
    
  }
  # define mean matrix for random effects
  mu[1] <- beta0
  mu[2] <- delta0
  mu[3] <- nu0
  
  beta0 ~ dnorm(a, b)
  delta0 ~ dnorm(c, d)
  nu0 ~ dnorm(e, f)
  
  # define precision matrix for random effects
  # Omega <- inverse(Sigma)
  
  Sigma[1, 1] <- pow(sigma.beta, 2)
  Sigma[2, 2] <- pow(sigma.delta, 2)
  Sigma[3, 3] <- pow(sigma.nu, 2)
  
  Sigma[1, 2] <- rho_bd * sigma.beta * sigma.delta
  Sigma[2, 1] <- Sigma[1, 2]
  Sigma[1, 3] <- rho_bn * sigma.beta * sigma.nu
  Sigma[3, 1] <- Sigma[1, 3]
  Sigma[2, 3] <- rho_dn * sigma.delta * sigma.nu
  Sigma[3, 2] <- Sigma[2, 3]
  
  # half-t prior on sigma.delta
  sigma.beta ~ dt(0, 1, 1) T(0,)
  sigma.nu ~ dt(0, 1, 1) T(0,)
  sigma.delta ~ dt(0, 1, 1) T(0,)
  
  # correlations
  rho_bd ~ dunif(-1, 1)
  rho_bn ~ dunif(-1, 1)
  rho_dn ~ dunif(-1, 1)
  
  

  #draw M new observations for each parameter w/ random effects
  
  for(j in 1:M){
    thetanew[j, 1:3] ~ dmnorm.vcov(mu[1:3], Sigma[1:3, 1:3])
    betanew[j] <- thetanew[j, 1]
    deltanew[j] <-  thetanew[j, 2]
    psinew[j] <- 1 / (1  + exp(-thetanew[j, 3]))
    
    pi11new[j] <- 1 / (1 + exp(-(betanew[j] + deltanew[j]/2))) * psinew[j]                  # P(event, risk)
    pi10new[j] <- (1 - 1 / (1 + exp(-(betanew[j] + deltanew[j] / 2)))) * psinew[j]          # P(no event, risk)
    pi01new[j] <- 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2))) * (1 - psinew[j])          # P(event, no risk)
    pi00new[j] <- (1 - 1 / (1 + exp(-(betanew[j] - deltanew[j] / 2)))) * (1 - psinew[j])    # P(no event, no risk)
    
    sensnew[j] <- pi11new[j] / (pi11new[j] + pi01new[j])
    specnew[j] <- pi00new[j] / (pi00new[j] + pi10new[j])
    
    PPVnew[j] <- pi11new[j] / (pi11new[j] + pi10new[j])
    NPVnew[j] <- pi00new[j] / (pi00new[j] + pi01new[j])
    
    LRpnew[j] <- sensnew[j] / (1 - specnew[j])
    LRmnew[j] <- (1 - sensnew[j]) / specnew[j]
    
  }
}", fill = TRUE)
sink()

### fun trial using multivariate normal instead of indep normals
S <- length(a) 
n.tot <- rowSums(n)
M <- 100

# inits

delta.init <- rep(0, S)
beta.init <- rep(-1, S)
nu.init <- rep(-1, S)

deltanew.init <- rep(0, M)
betanew.init <- rep(0, M)
nunew.init <- rep(0, M)

rho_bd.init <- 0.1
rho_bn.init <- 0
rho_dn.init <- -0.1


theta.init <- cbind(beta.init, delta.init, nu.init)
thetanew.init <- cbind(betanew.init, deltanew.init, nunew.init)

meta.data <- list(S = S, Y=Y, n=n, M = M, n.tot=n.tot, a=-2, b=0.5, c=1, d=0.5, 
                  e=0, f=0.5)



meta.multivar.inits <- rep(list(list(theta = theta.init, thetanew = thetanew.init,
                            beta0 = -1, nu0 = -1, delta0 = 0, 
                            sigma.delta = 0.5, sigma.beta = 0.5, sigma.nu = 0.5,
                            rho_bd = rho_bd.init, rho_bn = rho_bn.init, rho_dn = rho_dn.init)), 1)


meta.multivar.params <- c("rho_bd", "rho_bn", "rho_dn", "delta0", "nu0", "beta0")

meta.multivar <- jags(meta.data, meta.multivar.inits, meta.multivar.params, "meta_confusion_multivariate.txt", n.chains=1,
                   n.iter=21000, n.burnin=1000, n.thin=1, DIC=F)



# inits

delta.init <- rep(0, S)
beta.init <- rep(-1, S)
nu.init <- rep(-1, S)



deltanew.init <- rep(0, M)
betanew.init <- rep(0, M)
nunew.init <- rep(0, M)

meta.data <- list(S = S, Y=Y, n=n, M = M, n.tot=n.tot, a=-2, b=0.5, c=1, d=0.5, 
                  e=0, f=0.5)

init.gen <- function(){
  list(
    delta <- rnorm(S),
    beta <- rnorm(S),
    nu <- rnorm(S),
    delta0 <- rnorm(1),
    beta0 <- rnorm(1),
    nu0 <- rnorm(1),
    sigma.delta <- runif(1, 0.2, 1),
    sigma.beta <- runif(1, 0.2, 1),
    sigma.nu <- runif(1, 0.2, 1),
    deltanew <- rnorm(M),
    betanew <- rnorm(M),
    nunew <- rnorm(M)
  )
}


meta.params <- c("sensnew", "specnew", "PPVnew", "NPVnew", "LRpnew", "LRmnew")

meta_anal1 <- jags(meta.data, inits = init.gen, meta.params, "meta_confusion.txt", n.chains=3,
                  n.iter=21000, n.burnin=1000, n.thin=2, DIC=F)

meta_sims1 <- meta_anal1$BUGSoutput$sims.matrix
# 1-100 are pi00
# 101-200 are pi01
# 201-300 are pi10
# 301-400 are pi11

# sens1 <- apply(meta_sims1[,301:400] / (meta_sims1[,301:400] + meta_sims1[,101:200]), 
#               1, quantile, 0.5)
# 
# spec1 <- apply(meta_sims1[,1:100] / (meta_sims1[,1:100] + meta_sims1[,201:300]),
#               1, quantile, 0.5)
# 
# PPV1 <- apply(meta_sims1[,301:400] / (meta_sims1[,301:400] + meta_sims1[,201:300]),
#              1, quantile, 0.5)
# NPV1 <- apply(meta_sims1[,1:100] / (meta_sims1[,1:100] + meta_sims1[,101:200]),
#              1, quantile, 0.5)
# 
# LRp1 <- apply((meta_sims1[,301:400] * (meta_sims1[,1:100] + meta_sims1[,201:300])) / (meta_sims1[,201:300] * (meta_sims1[,301:400] + meta_sims1[,101:200])),
#              1, quantile, 0.5)
# 
# LRm1 <- apply((meta_sims1[,101:200] * (meta_sims1[,1:100] + meta_sims1[,201:300])) / (meta_sims1[,1:100] * (meta_sims1[,301:400] + meta_sims1[,101:200])),
#              1, quantile, 0.5)

LRm.med <- apply(meta_sims1[,1:100], 1, quantile, 0.5)
LRp.med <- apply(meta_sims1[,101:200], 1, quantile, 0.5)
NPV.med <- apply(meta_sims1[,201:300], 1, quantile, 0.5)
PPV.med <- apply(meta_sims1[,301:400], 1, quantile, 0.5)
sens.med <- apply(meta_sims1[,401:500], 1, quantile, 0.5)
spec.med <- apply(meta_sims1[,501:600], 1, quantile, 0.5)

baddies1 <- which(LRp1 > 100)


