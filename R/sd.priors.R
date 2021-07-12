######## Showing half-cauchy vs IG cdfs
library(extraDistr)

hc0.1 <- rhcauchy(10000, sigma = 0.1)
hc0.5 <- rhcauchy(10000, sigma = 0.5)
hc1.0 <- rhcauchy(10000, sigma = 1.0)
ig32 <- sqrt(rinvgamma(n = 10000, alpha = 3, beta = 2))
ig42 <- sqrt(rinvgamma(n = 10000, alpha = 4, beta = 2))

P1 <- ecdf(hc0.1)
P2 <- ecdf(hc0.5)
P3 <- ecdf(hc1.0)
P4 <- ecdf(ig32)
P5 <- ecdf(ig42)

P1(1.5)
P2(1.5)
P3(1.5)
P4(1.5)
P5(1.5)

plot(P1, xlim = c(0, 10))
lines(P2)
lines(P3)
lines(P4)
lines(P5)

