library(mvtnorm)
library(MCMCglmm)
library(ggplot2)
library(progress)

## ------------------------- Prior model for monovariate

## Parameters
K           <- 100000   # Number of iterations
nu          <- 1000     # nu parameter
alpha.mu    <- 0        # alpha.mu parameter
alpha.V     <- 1        # alpha.V parameter

## Prior for extended parameter Va
# Variance of eta
pr_var_eta  <- rIW(n = K, V = diag(1), nu = nu)[,1]
# Alpha
pr_alpha    <- rnorm(K, alpha.mu, sqrt(alpha.V))
# Resulting prior on Va
pr_var_a    <- (pr_alpha^2) * pr_var_eta

## Prior for residual variance
# Essentially useless for a binary trait, as fix = 1...
# but keeping it in case we want to remove fix = 1
pr_var_r <- rIW(n = K, V = diag(1), nu = 0.02, fix = 1)[,1]

## Resulting prior for heritability
pr_herit <- pr_var_a / (pr_var_a + pr_var_r)

## Plotting the density of the prior
qplot(x = pr_herit, geom = "density")

## ------------------------- Prior model for multivariate

## Parameters
K           <- 10000            # Number of iterations
nu          <- 1000             # nu parameter
dim         <- 2                # dimensions of the multivariate model
alpha.mu    <- c(0,0)           # alpha.mu parameter
alpha.V     <- diag(c(1, 1))    # alpha.V parameter

## Prior for extended parameter Va
# Variance of eta
pr_var_eta  <- rIW(n = K, V = diag(dim), nu = nu)
# Alpha
pr_alpha    <- rmvnorm(K, alpha.mu, alpha.V)
# Gather everything
pr_tot      <- cbind(pr_alpha, pr_var_eta)

## Resulting prior for G matrix
pr_var_G <- lapply(1:K, function(i){
    # Formatting alpha and eta matrices
    alpha <- diag(pr_alpha[i, ])
    ETA <- matrix(pr_var_eta[i, ], ncol = dim)
    # G is the following result of alpha and VCV of eta
    t(alpha) %*% ETA %*% alpha
})

## Prior for R matrix
# Note that the last dimension is binary as fix = dim
tmp <- rIW(n = K, V = diag(dim), nu = 2, fix = dim)
# Formatting a list
pr_var_R <- lapply(1:K, function(i){
    matrix(tmp[i, ], ncol = dim)
})

## Prior for heritabilities and correlations
# Setting up the output vectors
herit1  <- numeric(K)
herit2  <- numeric(K)
corG    <- numeric(K)
corE    <- numeric(K)
# Starting to loop over the samples
pb <- progress::progress_bar$new(total = K, force = TRUE, format  = "[:bar] :percent ETA: :eta")
for (i in 1:K){
  pb$tick()
  # Heritabilities
  herit1[i] <- pr_var_G[[i]][1, 1] / (pr_var_G[[i]][1, 1] + pr_var_R[[i]][1, 1])
  herit2[i] <- pr_var_G[[i]][2, 2] / (pr_var_G[[i]][2, 2] + pr_var_R[[i]][2, 2])
  # Genetic and env. correlations
  corG[i] <- cov2cor(pr_var_G[[i]])[1, 2]
  corE[i] <- cov2cor(pr_var_R[[i]])[1, 2]
}

## Plotting the densities
qplot(x = herit1, geom = "density")
qplot(x = herit2, geom = "density")
qplot(x = corG, geom = "density")
qplot(x = corE, geom = "density")
