## Figure 1, design example
setwd("/mypath")

library(R2OpenBUGS)
library(lattice)

## assume estimates were based on 14 urologic centers, with 50% serving low SES patients
ss = 14

n = ss
c = n / 2

p <- 0.50

## data
x <- c( rep( 1, p * n ), rep( 0, (1 - p) * n ) )

sum.xsq <- sum(x ^ 2)
sum.xminusxbarsq <- sum( (x - mean(x)) ^ 2 )

gam0.hat = -3.885
gam1.hat = 0.360
sigmasq.hat = 0.307 ^ 2
df = n - 2
sigmasq.u.hat = 0.104 ^ 2

## hypothesized difference to detect
delta = 0.15

## assume new trial with 10 pts per center, 20 centers (10 per arm)
m2 = 10
n2 = 20
c2 = n2 / 2

## r to openbugs
data_list <- list("gam0.hat", "gam1.hat", "sigmasq.hat", "sigmasq.u.hat", "df", "n", "c", "p", "delta",
                  "sum.xsq", "sum.xminusxbarsq", "m2", "c2")

inits <- function() {
    list( sigma = 0.5, sigma.u = 0.5, gam0 = -1.5, gam1 = 1 )
}

parameters <- c( "sigmasq", "sigmasq.u", "gam0", "gam1", "within", "power" )

niter <- 2000000
nthin <- 1
nburn <- 1900000
nchains <- 1

output = bugs(data = data_list, inits = inits, parameters = parameters,
              model.file = "obugs-fig1-designExample.txt", n.chains = nchains, n.thin = nthin,
              n.iter = niter, n.burnin = nburn, DIC = TRUE, codaPkg = TRUE,
              working.directory = "/mypath")

out.coda <- read.bugs(output)
summary(out.coda)
z <- as.matrix(out.coda)

out.density <- density( z[, 4] )
mode.posterior <- out.density$x[out.density$y == max(out.density$y)]
median <- median( z[, 4] )
mean <- mean( z[, 4] )
 
within <- ( 1 - p ) * exp( gam0.hat + sigmasq.hat / 2 ) + p * exp( gam0.hat + gam1.hat + sigmasq.hat / 2 )
freq.power <- pnorm( delta / ( sqrt( ( 2 / c2 ) * (within / m2 + sigmasq.u.hat) ) ) - 1.959964 )

results <- matrix(0, nrow = 1, ncol = 9)
results[1,] <- c( ss, c2, m2, quantile( z[, 4], 0.025 ), quantile( z[, 4], 0.10 ), mode.posterior, median, mean, freq.power )
colnames(results) <- c("N1", "c2", "m2", "power2.5", "power10.0", "power.mode", "power.median", "power.mean", "freq.power")

results

res <- data.frame(results)

# figure 1, posterior distribution power, 14 centers contributed to the estimates
power.dis <- z[, 4] 
hist( power.dis, prob = TRUE, breaks = 40, main = "Posterior Distribution for Power",
     xlab = "Power", col = "gray", lwd = 1 )
lines( density( power.dis ), col = "black", lwd = 2 )

out.density <- density( power.dis )
mode.posterior <- out.density$x[out.density$y == max(out.density$y)]
abline(v = median( power.dis ), lty = 2, col = 2)
abline(v = mean( power.dis ), lty = 2, col = "seagreen")
abline(v = quantile( power.dis, 0.10), lty = 2, col = 5)
abline(v = quantile( power.dis, 0.025), lty = 2, col = "blue")
legend("topleft", c("posterior median", "posterior mean", "10 percent quantile", "2.5 percent quantile"),
       col = c(2, "seagreen", 5, "blue"), lty = c(2, 2, 2, 2), lwd = c(1, 1, 1, 1),
       bty = "n", cex = 0.9 )

