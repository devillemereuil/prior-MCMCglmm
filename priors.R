library(mvtnorm)
library(MCMCglmm)
library(ggplot2)

##Writing the prior model for monovariate
K <- 100000	#Number of iterations
nu <- 1000	#nu parameter
alpha.mu <- 0	#alpha.mu parameter
alpha.V <- 1	#alpha.V parameter

#Prior for extended parameter Va
pr_var_eta <- rIW(n=K,V=diag(1),nu=nu)[,1]
pr_alpha <- rnorm(K,alpha.mu,sqrt(alpha.V))
pr_var_a <- (pr_alpha**2) * pr_var_eta

#Prior for residual variance
pr_var_r <- rIW(n=K,V=diag(1),nu=0.01)[,1]

#Resulting prior for heritability
pr_herit <- pr_var_a / (pr_var_a + 1)

#Plotting the density of the prior
qplot(x=pr_herit,geom="density")

##Writing the prior model for multi-response
K <- 100000		#Number of iterations
nu <- 1000		#nu parameter
dim <- 2		#dimensions of the multivariate model
alpha.mu <- c(0,0)	#alpha.mu parameter
alpha.V <- diag(dim)	#alpha.V parameter

#Prior for extended parameter Va
pr_var_eta <- rIW(n=K,V=diag(dim),nu=nu)
pr_alpha <- rmvnorm(K,alpha.mu,alpha.V)
pr_tot <- cbind(pr_alpha,pr_var_eta)

#Prior for G matrix
pr_var_G <- lapply(1:K,function(i){
  vec <- pr_tot[i,]
  alpha <- diag(c(vec[1],vec[2]))
  ETA <- matrix(c(vec[3],vec[4],vec[5],vec[6]),ncol=dim)
  t(alpha) %*% ETA %*% alpha
})

#Prior for R matrix (second dimension is binary)
tmp <- rIW(n=K,V=diag(dim),nu=0.1,fix=2)
pr_var_R <- lapply(1:K,function(i){
  matrix(tmp[i,],ncol=dim)
})

#Prior for heritabilities and correlations
herit1 <- numeric(K)
herit2 <- numeric(K)
corG <- numeric(K)
corE <- numeric(K)
for (i in 1:K){
  print(i)
  herit1[i] <- pr_var_G[[i]][1,1]/(pr_var_G[[i]][1,1] + pr_var_R[[i]][1,1])
  herit2[i] <- pr_var_G[[i]][2,2]/(pr_var_G[[i]][2,2] + pr_var_R[[i]][2,2])
  corG[i] <- cov2cor(pr_var_G[[i]])[1,2]
  corE[i] <- cov2cor(pr_var_R[[i]])[1,2]
}

#Plotting the densities
qplot(x=herit1,geom="density")
qplot(x=herit2,geom="density")
qplot(x=corG,geom="density")
qplot(x=corE,geom="density")