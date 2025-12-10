library(tidyverse)
library(tidytuesdayR)

# Load and prepare data
rm(list = ls())
set.seed(432)
gas_data <- tt_load(2025, week = 26)
gas_data <- gas_data$weekly_gas_prices
regular <- gas_data |>
  dplyr::filter(fuel == "gasoline", grade == "all", formulation == "all") |>
  dplyr::select(date, price)

# Prepare design matrix - ONLY intercept and time
T <- nrow(regular)
t.scaled <- (1:T) / T
X <- cbind(1, t.scaled)
y <- regular$price

hts.mcmc <- function(y, X, priors, inits, n.iter){
  
  #
  #  Fits Hierarchical Time Series Model
  #  y_t ~ N(x_t'beta + z_t, sigma^2_y)
  #  z_t ~ N(rho*z_{t-1}, sigma^2_z)
  #
  
  ####
  ####  Setup Variables
  ####
  
  T <- length(y)
  p <- ncol(X)
  
  ####
  ####  Priors and Starting values
  ####
  
  mu.beta <- priors$mu.beta
  Sigma.beta <- priors$Sigma.beta
  Sigma.beta.inv <- solve(Sigma.beta)
  beta.est <- inits$beta
  
  q.y <- priors$q.y
  r.y <- priors$r.y
  s2.y.est <- inits$s2.y
  
  q.z <- priors$q.z
  r.z <- priors$r.z
  s2.z.est <- inits$s2.z
  
  a.rho <- priors$a.rho
  b.rho <- priors$b.rho
  rho.est <- inits$rho
  
  mu.0 <- priors$mu.0
  s2.0 <- priors$s2.0
  z.est <- inits$z
  
  # M-H tuning
  rho.tune <- 0.2
  rho.accept <- 0
  
  n.burn <- 0.4*n.iter
  thin <- 10
  idx.save <- seq(n.burn, n.iter, by=thin)
  beta.save <- matrix(NA, length(idx.save), p)
  s2.y.save <- rep(NA, length(idx.save))
  s2.z.save <- rep(NA, length(idx.save))
  rho.save <- rep(NA, length(idx.save))
  z.save <- matrix(NA, length(idx.save), T)
  
  ####
  ####  Begin MCMC Loop
  ####
  
  for(k in 1:n.iter){
    if(k%%1000==0) cat(k, " ")
    
    ####
    ####  Sample beta
    ####
    
    A <- t(X) %*% X / s2.y.est + Sigma.beta.inv
    b <- t(X) %*% (y - z.est) / s2.y.est + Sigma.beta.inv %*% mu.beta
    A.inv <- solve(A)
    mu.tmp <- A.inv %*% b
    beta.est <- c(mu.tmp + t(chol(A.inv)) %*% rnorm(p))
    
    ####
    ####  Sample s2.y
    ####
    
    q.tmp <- T / 2 + q.y
    resid.y <- y - X %*% beta.est - z.est
    r.tmp <- sum(resid.y^2) / 2 + 1 / r.y
    s2.y.est <- 1 / rgamma(1, q.tmp, r.tmp)
    
    ####
    ####  Sample s2.z
    ####
    
    q.tmp <- T / 2 + q.z
    resid.z <- c(z.est[1] - mu.0, z.est[2:T] - rho.est * z.est[1:(T-1)])
    r.tmp <- sum(resid.z^2) / 2 + 1 / r.z
    s2.z.est <- 1 / rgamma(1, q.tmp, r.tmp)
    
    ####
    ####  Sample rho (Metropolis-Hastings)
    ####
    
    logit.rho <- log(rho.est / (1 - rho.est))
    logit.rho.star <- rnorm(1, logit.rho, rho.tune)
    rho.star <- exp(logit.rho.star) / (1 + exp(logit.rho.star))
    
    ll.cur <- dnorm(z.est[1], mu.0, sqrt(s2.0), log = TRUE) + 
      sum(dnorm(z.est[2:T], rho.est * z.est[1:(T-1)], sqrt(s2.z.est), log = TRUE))
    ll.star <- dnorm(z.est[1], mu.0, sqrt(s2.0), log = TRUE) + 
      sum(dnorm(z.est[2:T], rho.star * z.est[1:(T-1)], sqrt(s2.z.est), log = TRUE))
    
    lp.cur <- dbeta(rho.est, a.rho, b.rho, log = TRUE)
    lp.star <- dbeta(rho.star, a.rho, b.rho, log = TRUE)
    
    lj.cur <- log(rho.est) + log(1 - rho.est)
    lj.star <- log(rho.star) + log(1 - rho.star)
    
    log.mh <- (ll.star + lp.star + lj.star) - (ll.cur + lp.cur + lj.cur)
    
    if(log(runif(1)) < log.mh){
      rho.est <- rho.star
      if(k > n.burn) rho.accept <- rho.accept + 1
    }
    
    ####
    ####  Sample z
    ####
    
    # Sample z_0 (initial state)

    A.tmp <- rho.est^2 / s2.z.est + 1 / s2.0
    b.tmp <- rho.est * z.est[1] / s2.z.est + mu.0 / s2.0
    s2.tmp <- 1 / A.tmp
    mu.tmp <- s2.tmp * b.tmp
    z.est[1] <- rnorm(1, mu.tmp, sqrt(s2.tmp))
    
    # Sample z_t for t = 2, ..., T-1 (middle states)
    
    for(t in 2:(T-1)){
      A.tmp <- 1 / s2.y.est + (rho.est^2 + 1) / s2.z.est
      b.tmp <- (y[t] - sum(X[t,] * beta.est)) / s2.y.est + 
        rho.est * (z.est[t+1] + z.est[t-1]) / s2.z.est
      s2.tmp <- 1 / A.tmp
      mu.tmp <- s2.tmp * b.tmp
      z.est[t] <- rnorm(1, mu.tmp, sqrt(s2.tmp))
    }
    
    # Sample z_T (final state)
    
    A.tmp <- 1 / s2.y.est + 1 / s2.z.est
    b.tmp <- (y[T] - sum(X[T,] * beta.est)) / s2.y.est + 
      rho.est * z.est[T-1] / s2.z.est
    s2.tmp <- 1 / A.tmp
    mu.tmp <- s2.tmp * b.tmp
    z.est[T] <- rnorm(1, mu.tmp, sqrt(s2.tmp))
    
    ####
    ####  Save Samples
    ####
    
    if(k %in% idx.save){
      beta.save[which(idx.save==k), ] <- beta.est
      s2.y.save[which(idx.save==k)] <- s2.y.est
      s2.z.save[which(idx.save==k)] <- s2.z.est
      rho.save[which(idx.save==k)] <- rho.est
      z.save[which(idx.save==k), ] <- z.est
    }
  }
  cat("\n")
  
  accept.rate <- rho.accept / length((n.burn+1):n.iter)
  cat("Rho acceptance rate:", round(accept.rate, 3), "\n")
  
  ####
  ####  Write Output
  ####
  
  list('beta'=beta.save,
       's2.y'=s2.y.save,
       's2.z'=s2.z.save,
       'rho'=rho.save,
       'z'=z.save,
       'accept.rate'=accept.rate)
}

####
#### Set Priors
####

priors <- list(
  mu.beta = c(0, 0),
  Sigma.beta = diag(100, 2),
  q.y = 2,
  r.y = 0.01,
  q.z = 2,
  r.z = 0.01,
  a.rho = 20,
  b.rho = 2,
  mu.0 = 0,
  s2.0 = 1
)

####
#### Set Initial Values
####

inits <- list(
  beta = c(mean(y), 0),
  s2.y = var(y) * 0.5,
  s2.z = var(y) * 0.5,
  rho = 0.9,
  z = rep(0, T)
)

####
#### Run MCMC
####

set.seed(123)
out <- hts.mcmc(y, X, priors, inits, n.iter = 20000)

####
#### Posterior Summary Statistics
####

# Create summary table
params <- c("Intercept", "Time Trend", "Sigma^2_y", "Sigma^2_z", "Rho")
samples <- list(out$beta[,1], out$beta[,2], out$s2.y, out$s2.z, out$rho)

summary_table <- data.frame(
  Parameter = params,
  Mean = sapply(samples, mean),
  SD = sapply(samples, sd),
  Q2.5 = sapply(samples, quantile, 0.025),
  Q50 = sapply(samples, quantile, 0.50),
  Q97.5 = sapply(samples, quantile, 0.975)
)

rownames(summary_table) <- NULL

cat("========================================\n")
cat("   POSTERIOR PARAMETER ESTIMATES\n")
cat("========================================\n\n")
print(summary_table, digits=4, row.names=FALSE)
cat("\n")

####
#### Trace Plots
####

par(mfrow=c(2,3), mar=c(4,4,2,1))

plot(out$beta[,1], type='l', main='Trace: Intercept', 
     xlab='Iteration', ylab='Beta_0')
abline(h=mean(out$beta[,1]), col='red', lwd=2)

plot(out$beta[,2], type='l', main='Trace: Time Trend', 
     xlab='Iteration', ylab='Beta_1')
abline(h=mean(out$beta[,2]), col='red', lwd=2)

plot(out$s2.y, type='l', main='Trace: Sigma^2_y', 
     xlab='Iteration', ylab='Sigma^2_y')
abline(h=mean(out$s2.y), col='red', lwd=2)

plot(out$s2.z, type='l', main='Trace: Sigma^2_z', 
     xlab='Iteration', ylab='Sigma^2_z')
abline(h=mean(out$s2.z), col='red', lwd=2)

plot(out$rho, type='l', main='Trace: Rho', 
     xlab='Iteration', ylab='Rho')
abline(h=mean(out$rho), col='red', lwd=2)

####
#### Fitted vs Observed
####

par(mfrow=c(1,1), mar=c(4,4,3,1))
z.mean <- apply(out$z, 2, mean)
beta.mean <- apply(out$beta, 2, mean)
fitted.mean <- X %*% beta.mean + z.mean

plot(regular$date, y, type='l', col='black', lwd=1,
     xlab='Date', ylab='Gas Price ($)',
     main='Observed vs Fitted Gas Prices')
lines(regular$date, fitted.mean, col='red', lwd=2)
legend('topleft', c('Observed', 'Fitted'), 
       col=c('black', 'red'), lty=1, lwd=c(1,2), bty='n')