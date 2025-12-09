library(tidyverse)
library(tidytuesdayR)

rm(list = ls())

set.seed(432)

gas_data <- tt_load(2025, week = 26)
gas_data <- gas_data$weekly_gas_prices

regular <- gas_data |>
  filter(fuel == "gasoline", grade == "all", formulation == "all") |>
  select(date, price)

y <- regular$price
X <- matrix(1, nrow = length(regular$date), ncol = 2)
X[,2] <- regular$date

priors <- list(
  "mu.beta" <- coef(lm(price ~ date, data = regular)),
  "Sig.beta" <- 1E3 * diag(2),
  "q.y" <- 1E-3,
  "r.y" <- 1E3,
  "q.z" <- 1E-3,
  "r.z" <- 1E3,
  "a" <- 1,
  "b" <- 1
)

inits <- list(
  "beta" <- coef(lm(price ~ date, data = regular)),
  "s2.y" <- var(y),
  "rho" <- 0,
  "s2.z" <- 10,
  "mu.0" <- 0,
  "s2.0" <- 1E3
)

n.iter <- 1E4

time.series.mcmc <- function(y, X, priors, inits, n.iter) {
  mu.beta <- priors$mu.beta
  Sig.beta <- priors$Sig.beta
  beta <- inits$beta
  
  q.y <- priors$q.y
  r.y <- priors$r.y
  s2.y <- inits$s2.y
  
  mu.0 <- priors$mu.0
  s2.0 <- priors$s2.0
  z.0 <- rnorm(1, mu.0, s2.0)
  
  q.z <- priors$q.z
  r.z <- priors$r.z
  s2.z <- inits$s2.z
  
  a <- priors$a
  b <- priors$b
  rho <- inits$rho
  z.t_1 <- z.0
  z.t <- rnorm(1, rho * z.t_1, s2.z)
}

out <- time.series.mcmc(y, X, priors, inits, n.iter)
