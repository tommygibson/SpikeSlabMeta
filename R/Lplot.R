########## Show how changing L (# of draws from predictive dist'n in each iteration of MCMC)
########## changes inferences


library(R2jags)
library(here)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(xtable)

# font_import(pattern = "lmroman*")

source("R/spike.functions.R")

#### Make plot for a given set of hyperparameters

delta0 <- 2
beta0 <- log(.15 / .85)
nu0 <- log(.15 / .85)

sigma.d <- sigma.b <- sigma.n <- c(0.5, 1)

S <- 10
N <- 500
L <- c(1, 10, 100, 1000)

LRp <- list()
for(i in 1:2) LRp[[i]] <- matrix(nrow = 10000, ncol = length(L))


for(j in 1:length(sigma.d)){
  
  set.seed(1516)
  
  delta.i <- rnorm(S, mean = delta0, sd = sigma.d[j])
  beta.i <- rnorm(S, mean = beta0, sd = sigma.b[j])
  nu.i <- rnorm(S, mean = nu0, sd = sigma.n[j])
  
  pi11.i <- 1 / (1 + exp(-(beta.i + delta.i / 2))) * 1 / (1 + exp(-nu.i))
  pi10.i <- (1 - 1 / (1 + exp(-(beta.i + delta.i / 2)))) * 1 / (1 + exp(-nu.i))
  pi01.i <- 1 / (1 + exp(-(beta.i - delta.i / 2))) * (1 - 1 / (1 + exp(-nu.i)))
  pi00.i <- (1 - 1 / (1 + exp(-(beta.i - delta.i / 2)))) * (1 - 1 / (1 + exp(-nu.i)))
  
  probs <- cbind(pi11.i, pi10.i, pi01.i, pi00.i) # collect all the cell probabilities
  
  tables <- t(apply(probs, 1, rmultinom, n = 1, size = N)) # generate contingency tables
  
  y <- tables[,c(1, 3)]
  n <- cbind(rowSums(tables[, 1:2]), rowSums(tables[, 3:4]))
  n.tot <- rep(N, S)
  
  for(i in 1:length(L)){
    
    # prepare data for going into jags
    
    M <- L[i]
    
    
    meta.dat <- list(M = M, S = S, n = n, y = y, n.tot = n.tot,
                     a = -1, b = 0.5, c = 0, d = 0.5, e = -1, f = 0.5)
    
    meta.params <- c("LRpnew")
    
    fit <- jags(data = meta.dat, inits = init.gen, parameters.to.save = meta.params,
                model.file = "R/meta_confusion.txt",
                n.chains = 2, n.iter = 11000, n.thin = 2, n.burnin = 1000, DIC = FALSE)
    
    if(M < 2) {
      LRp[[j]][, i] <- fit$BUGSoutput$sims.matrix
    }
    else {
      LRp[[j]][, i] <- apply(fit$BUGSoutput$sims.matrix, 1, mean)
    }
    
    
  }
}


LRp.low <- cbind.data.frame(as.vector(LRp[[1]]),
                            as.factor(as.character(rep(L, each = dim(LRp[[1]])[1]))))
LRp.high <- cbind.data.frame(as.vector(LRp[[2]]),
                            as.factor(as.character(rep(L, each = dim(LRp[[2]])[1]))))

names(LRp.low) <- names(LRp.high) <- c("LR+", "L")


LRp.lowhigh <- cbind.data.frame(rbind(LRp.low, LRp.high), 
                                as.factor(rep(c("SDs = 0.5", "SDs = 1"), each = dim(LRp.low)[1])))
names(LRp.lowhigh) <- c("LR+", "L", "SDs")

Lplot <- LRp.lowhigh %>%
           ggplot(mapping = aes(x = `LR+`, linetype = L, size = L)) +
             theme_bw() +
             theme(text = element_text(size = 12, family = "LM Roman 10")) +
             geom_density(fill = "black", alpha = 0.05) +
             facet_wrap( ~ SDs) +
             xlim(c(0, 10)) +
             scale_linetype_manual(values = c("dotted", "dotdash", "solid", "dashed")) +
             scale_size_manual(values = c(0.5, 0.5, 0.5, 0.5)) + 
             labs(ylab = c("", ""))


ggsave("TeX/L_plot.pdf", plot = Lplot, device = "pdf",
       width = 5, height = 5, units = "in")

simsum <- function(x){
  sd <- apply(x, 2, sd)
  length <- as.vector(diff(apply(x, 2, quantile, c(.025, .975))))
  
  return(sigfig(cbind(sd, length), n = 3))
}

sd0.5 <- simsum(LRp[[1]])
sd1 <- simsum(LRp[[2]])

Lplot.tab <- cbind.data.frame(sd0.5, sd1)
rownames(Lplot.tab) <- c("L = 1", "L = 10", "L = 100", "L = 1000")

print(xtable(Lplot.tab, caption = "Posterior distribution of LR+ for varying L for two simulated datasets",
             type = "latex"),
      file = "TeX/Lplot.tex")


                            