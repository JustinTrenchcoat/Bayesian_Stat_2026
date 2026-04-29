# this is the testing page
# in here we would only use data to fit global OLS
# and use them to set up groups
# we will use the global OLS for our Bayesian Prior



# add OLS beta vs group size
# effective sample size/MC error vs mean of THETA or other ones
# beta j ols - beta j posterior vs sample size
data <- read.csv("data/grouped.csv")

# -----------------------------
# Data preprocessing
# -----------------------------
data$fuel_type <- factor(data$fuel_type)
data$brand     <- factor(data$brand)
data$accident  <- factor(data$accident)

split_ids <- split(seq_len(nrow(data)), data$brand)

train_idx <- unlist(
  lapply(split_ids, function(idx){
    
    n <- length(idx)
    
    if(n <= 2){
      return(idx)
    }
    
    sample(idx, floor(0.8 * n))
  })
)

data_train <- data[train_idx, ]
data_test  <- data[-train_idx, ]

#standardizing the numerical variables
m_mean <- mean(data_train$mileage)
m_sd   <- sd(data_train$mileage)

a_mean <- mean(data_train$age)
a_sd   <- sd(data_train$age)

data_train$mileage_z <- (data_train$mileage - m_mean) / m_sd
data_test$mileage_z  <- (data_test$mileage - m_mean) / m_sd

data_train$age_z <- (data_train$age - a_mean) / a_sd
data_test$age_z  <- (data_test$age - a_mean) / a_sd

# -----------------------------
# Global pooled OLS
# -----------------------------
fit <- lm(price ~ mileage_z + fuel_type + age_z + accident, data = data_train)
summary(fit)

# -----------------------------
# Group setup
# -----------------------------
ids <- levels(data_train$brand)
m   <- length(ids)

X_full <- model.matrix(~ mileage_z + fuel_type + age_z + accident, data = data_train)
Y_full <- data_train$price

Y <- vector("list", m)
X <- vector("list", m)
N <- numeric(m)

for(j in 1:m){
  idx <- which(data_train$brand == ids[j])
  
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
    
    if(s %% 1000 == 0) cat(s, s2, "\n")
    
    S2.b    <- c(S2.b, s2)
    THETA.b <- rbind(THETA.b, as.numeric(theta))
    
    Sigma_now <- solve(iSigma)
    
    Sigma.ps <- Sigma.ps + Sigma_now
    BETA.ps  <- BETA.ps + BETA
    
    SIGMA.PS <- rbind(SIGMA.PS, c(Sigma_now))
    BETA.pp  <- rbind(BETA.pp, rmvnorm(1, c(theta), Sigma_now))
    BETA.store[[length(BETA.store)+1]] <- BETA
  }
}

# posterior means
Sigma.ps <- Sigma.ps / length(S2.b)
BETA.ps  <- BETA.ps / length(S2.b)

# -----------------------------
# Diagnostics
# -----------------------------
library(coda)

# effective size checks
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(BETA.pp[,1])

# acf checks
acf(S2.b)
acf(THETA.b[,1])
acf(BETA.pp[,1])

# Visualization of posterior coefficients
# customize brand, BMW is subject to change to any brand in the data set:
j <- which(ids=="Acura")
# customize coefficient search, two variables are subject to change:
k1 <- which(colnames(BETA.store[[1]]) == "mileage_z")
k2 <- which(colnames(BETA.store[[1]]) == "age_z")
# extract data from stored coefficients
beta_x <- sapply(BETA.store, function(mat) mat[j, k1])
beta_y <- sapply(BETA.store, function(mat) mat[j, k2])
# trace plot, title and variable subject to change
plot(beta_x, type="l",
     main="Coefficient Trace Plot",
     ylab="beta", xlab="Iteration")
# joint posterior distribution
plot(beta_x, beta_y,
     pch=19, cex=.5,
     xlab="[subject to change] Coefficient",
     ylab="[subject to change] Coefficient",
     main="[brand name] Posterior Joint Draws")

######################################################
# prediction performance assessment:
RSE <- list()
test_length <- nrow(data_test)
for (d in 1:test_length){
  
}
for (d in 1:test_length){
  idx <- which(data_test$brand[d])
  beta39 <- BETA.ps[39, ]
  predictions <- X_test[idx, ] %*% beta39
  
}