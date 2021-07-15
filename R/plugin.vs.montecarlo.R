nui <- rnorm(n = 1000, -2, 0.75)
betai <- rnorm(n = 1000, -2, 0.75)
deltai <- rnorm(n = 1000, 2, 0.75)
nu0 <- mean(nui)
beta0 <- mean(betai)
delta0 <- mean(deltai)

psi <- 1 / (1 + exp(-nui))
pi1 <- 1 / (1 + exp(-(betai + deltai / 2)))
pi0 <- 1 / (1 + exp(-(betai - deltai / 2)))

psi.plug <- 1 / (1 + exp(-nu0))
pi1.plug <- 1 / (1 + exp(-(beta0 + delta0 / 2)))
pi0.plug <- 1 / (1 + exp(-(beta0 - delta0 / 2)))

sens <- pi1 * psi / (pi1 * psi + pi0 * (1 - psi))
spec <- (1 - pi0) * (1 - psi) / ((1 - pi0) * (1 - psi) + (1 - pi1) * psi)

sens.plug <- pi1.plug * psi.plug / (pi1.plug * psi.plug + pi0.plug * (1 - psi.plug))
spec.plug <- (1 - pi0.plug) * (1 - psi.plug) / ((1 - pi0.plug) * (1 - psi.plug) + (1 - pi1.plug) * psi.plug)

LRP <- sens / (1 - spec)
LRM <- (1 - sens) / spec

LRP.plug <- sens.plug / (1 - spec.plug)
LRM.plug <- (1 - sens.plug) / spec.plug

cbind(apply(cbind(psi, pi1, pi0, sens, spec, LRP, LRM), 2, mean), c(psi.plug, pi1.plug, pi0.plug, sens.plug, spec.plug, LRP.plug, LRM.plug))

########## 

