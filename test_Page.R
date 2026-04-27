data <- read.csv("data/cleaned.csv")

# -----------------------------
# Data preprocessing
# -----------------------------
data$fuel_type <- factor(data$fuel_type)
data$brand     <- factor(data$brand)
data$accident  <- factor(data$accident)

data$mileage <- as.numeric(scale(data$mileage))
data$age     <- as.numeric(scale(data$age))

# -----------------------------
# Global pooled OLS
# -----------------------------
fit <- lm(price ~ mileage + fuel_type + age + accident, data = data)
summary(fit)

# -----------------------------
# Group setup
# -----------------------------
ids <- levels(data$brand)
m   <- length(ids)

X_full <- model.matrix(~ mileage + fuel_type + age + accident, data = data)
Y_full <- data$price

Y <- vector("list", m)
X <- vector("list", m)
N <- numeric(m)

for(j in 1:m){
  idx <- which(data$brand == ids[j])
  
  Y[[j]] <- Y_full[idx]
  X[[j]] <- X_full[idx, , drop = FALSE]
  N[j]   <- length(idx)
}

# -----------------------------
# Utilities
# -----------------------------
rmvnorm <- function(n, mu, Sigma){
  p <- length(mu)
  Z <- matrix(rnorm(n*p), n, p)
  Z %*% chol(Sigma) + matrix(mu, n, p, byrow = TRUE)
}

rwish <- function(n, nu0, S0){
  p <- nrow(S0)
  C <- chol(S0)
  out <- array(0, dim = c(p,p,n))
  
  for(i in 1:n){
    Z <- matrix(rnorm(nu0*p), nu0, p) %*% C
    out[,,i] <- t(Z) %*% Z
  }
  
  out
}

# -----------------------------
# Bayesian setup
# -----------------------------
beta_hat   <- coef(fit)
V_beta     <- vcov(fit)
sigma_hat2 <- summary(fit)$sigma^2

p <- length(beta_hat)

theta <- mu0 <- as.matrix(beta_hat, ncol = 1)

BETA <- matrix(rep(beta_hat, each = m), nrow = m, ncol = p, byrow = TRUE)
colnames(BETA) <- names(beta_hat)

s2  <- s20 <- sigma_hat2

nu0  <- 1
eta0 <- p + 2

Sigma <- S0 <- L0 <- V_beta*10
iL0    <- solve(L0)
iSigma <- solve(Sigma)

# -----------------------------
# Storage
# -----------------------------
THETA.b  <- NULL
S2.b     <- NULL
SIGMA.PS <- NULL
BETA.pp  <- NULL

Sigma.ps <- matrix(0, p, p)
BETA.ps  <- matrix(0, m, p)

BETA.store <- list()

# -----------------------------
# MCMC
# -----------------------------
set.seed(632)

nsim   <- 50000
thin   <- 10

for(s in 1:nsim){
  
  # update beta_j
  for(j in 1:m){
    
    XtX <- t(X[[j]]) %*% X[[j]]
    XtY <- t(X[[j]]) %*% Y[[j]]
    
    Vj <- solve(iSigma + XtX / s2)
    Ej <- Vj %*% (iSigma %*% theta + XtY / s2)
    
    BETA[j, ] <- rmvnorm(1, c(Ej), Vj)
  }
  
  # update theta
  Lm  <- solve(iL0 + m * iSigma)
  mum <- Lm %*% (iL0 %*% mu0 + iSigma %*% matrix(colSums(BETA), ncol = 1))
  
  theta <- t(rmvnorm(1, c(mum), Lm))
  
  # update Sigma
  mtheta <- matrix(as.numeric(theta), m, p, byrow = TRUE)
  
  SS <- S0 + t(BETA - mtheta) %*% (BETA - mtheta)
  
  iSigma <- rwish(1, eta0 + m, solve(SS))[,,1]
  
  # update s2
  RSS <- 0
  for(j in 1:m){
    resid_j <- Y[[j]] - X[[j]] %*% BETA[j, ]
    RSS <- RSS + sum(resid_j^2)
  }
  
  s2 <- 1 / rgamma(
    1,
    shape = (nu0 + sum(N))/2,
    rate  = (nu0*s20 + RSS)/2
  )
  
  # store
  if(s %% thin == 0){
    
    BETA.store[[length(BETA.store)+1]] <- BETA
    
    if(s %% 1000 == 0) cat(s, s2, "\n")
    
    S2.b    <- c(S2.b, s2)
    THETA.b <- rbind(THETA.b, as.numeric(theta))
    
    Sigma_now <- solve(iSigma)
    
    Sigma.ps <- Sigma.ps + Sigma_now
    BETA.ps  <- BETA.ps + BETA
    
    SIGMA.PS <- rbind(SIGMA.PS, c(Sigma_now))
    BETA.pp  <- rbind(BETA.pp, rmvnorm(1, c(theta), Sigma_now))
  }
}

# posterior means
Sigma.ps <- Sigma.ps / length(S2.b)
BETA.ps  <- BETA.ps / length(S2.b)

# -----------------------------
# Diagnostics
# -----------------------------
library(coda)

effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(BETA.pp[,1])

acf(S2.b)
acf(THETA.b[,1])
