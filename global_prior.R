# this is the page where the priors are set up from a global OLS model.
#-----------------------------
# Necessary library
#-----------------------------
library(coda)
# -----------------------------
# Data pre-processing
# -----------------------------
data <- read.csv("data/grouped.csv")
data$fuel_type <- factor(data$fuel_type)
data$brand     <- factor(data$brand)
data$accident  <- factor(data$accident)

#standardizing the numerical variables
# mileage
m_mean <- mean(data$mileage)
m_sd   <- sd(data$mileage)
data$mileage <- (data$mileage - m_mean) / m_sd

# age
a_mean <- mean(data$age)
a_sd   <- sd(data$age)
data$age <- (data$age - a_mean) / a_sd
# pairs plot for checking correlation, but not included in the report
vars <- data[, c("price", "mileage", "age", "accident")]
pairs(vars)
# -----------------------------
# Global pooled OLS
# -----------------------------
fit <- lm(price ~ mileage+ fuel_type + age + accident, data = data)
summary(fit)
res <- residuals(fit)
qqnorm(res, main = "Q-Q Plot of Residuals")
qqline(res, col = "red", lwd = 2)
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
# Functions 
# -----------------------------
rmvnorm <- function(n, mu, Sigma){
  mu_l <- length(mu)
  E <- matrix(rnorm(n*mu_l),n,mu_l)
  t(t(E%*% chol(Sigma))+c(mu))         
}

rwish <- function(n, nu0, S0){
  sS0 <- chol(S0)
  S <- array(dim=c(dim(S0),n))
  for (i in 1:n){
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1] ) %*% sS0
    S[,,i] <- t(Z)%*%Z
  }
  S[,,1:n]
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
BETA.store  <- NULL

Sigma.ps <- matrix(0, p, p)
BETA.ps  <- matrix(0, m, p)

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
    
    BETA[j, ] <- rmvnorm(1, Ej, Vj)
  }
  
  # update theta
  Lm  <- solve(iL0 + m * iSigma)
  mum <- Lm %*% (iL0 %*% mu0 + iSigma %*% matrix(colSums(BETA), ncol = 1))
  theta <- t(rmvnorm(1, mum, Lm))
  
  # update Sigma
  mtheta <- matrix(theta, m, p, byrow = TRUE)
  SS <- S0 + t(BETA - mtheta) %*% (BETA - mtheta)
  iSigma <- rwish(1, eta0 + m, solve(SS))
  
  # update s2
  RSS <- 0
  for(j in 1:m){
    resid_j <- Y[[j]] - X[[j]] %*% BETA[j, ]
    RSS <- RSS + sum(resid_j^2)
  }
  
  s2 <- 1 / rgamma(
    1,
    (nu0 + sum(N))/2,
    (nu0*s20 + RSS)/2
  )
  
  # store
  # # testing the thinning and ESS relation
  # if(s %% 1000 == 0) cat(s, s2, "\n")
  # S2.b    <- c(S2.b, s2)
  # THETA.b <- rbind(THETA.b, t(theta))
  # Sigma_now <- solve(iSigma)
  # Sigma.ps <- Sigma.ps + Sigma_now
  # BETA.ps  <- BETA.ps + BETA
  # SIGMA.PS <- rbind(SIGMA.PS, c(Sigma_now))

  if(s %% thin == 0){

    if(s %% 1000 == 0) cat(s, s2, "\n")
    S2.b    <- c(S2.b, s2)
    THETA.b <- rbind(THETA.b, t(theta))
    Sigma_now <- solve(iSigma)
    Sigma.ps <- Sigma.ps + Sigma_now
    BETA.ps  <- BETA.ps + BETA
    SIGMA.PS <- rbind(SIGMA.PS, c(Sigma_now))
    BETA.store <- rbind(BETA.store, BETA)
    }
}

# posterior means
Sigma.ps <- Sigma.ps / length(S2.b)
BETA.ps  <- BETA.ps / length(S2.b)

# Diagnostics
# -----------------------------
# effective size checks
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(THETA.b[,2])
effectiveSize(THETA.b[,3])
effectiveSize(THETA.b[,4])
effectiveSize(THETA.b[,5])
effectiveSize(THETA.b[,6])
effectiveSize(THETA.b[,7])
effectiveSize(THETA.b[,8])

# acf checks
acf(S2.b)
acf(THETA.b[,1])
acf(THETA.b[,2])
acf(THETA.b[,3])
acf(THETA.b[,4])
acf(THETA.b[,5])
acf(THETA.b[,6])
acf(THETA.b[,7])
acf(THETA.b[,8])

# trace plots 
# plot(S2.b, type = "l",
#      xlab = "Iteration",
#      ylab = expression(sigma^2),
#      main = "Trace plot for sigma^2")
# abline(h = mean(S2.b), col = "red", lwd = 2)
# plot(THETA.b[,1], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta 1")
# abline(h = mean(THETA.b[,1]), col = "red", lwd = 2)
# plot(THETA.b[,2], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta 2")
# abline(h = mean(THETA.b[,2]), col = "red", lwd = 2)
# plot(THETA.b[,3], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta in E85 Flex Fuel")
# abline(h = mean(THETA.b[,3]), col = "red", lwd = 2)
# plot(THETA.b[,4], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta 4")
# abline(h = mean(THETA.b[,4]), col = "red", lwd = 2)
# plot(THETA.b[,6], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta 6")
# abline(h = mean(THETA.b[,6]), col = "red", lwd = 2)
# plot(THETA.b[,7], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta 7")
# abline(h = mean(THETA.b[,7]), col = "red", lwd = 2)
# plot(THETA.b[,8], type = "l",
#      xlab = "Iteration",
#      ylab = expression(theta),
#      main = "Trace plot for Theta 8")
# abline(h = mean(THETA.b[,8]), col = "red", lwd = 2)
# Plot example in the report
plot(THETA.b[,5], type = "l",
     xlab = "Iteration",
     ylab = expression(theta),
     main = "Trace plot for Theta Hybrid")
abline(h = mean(THETA.b[,5]), col = "red", lwd = 2)
# -----------------------------

# Visualizations
# -----------------------------
# density comparison
# plot(density(THETA.b[,2],adj=2),xlim=range(BETA.store[,2]), 
#      main="",xlab="beta on mileage",ylab="posterior density",lwd=2)
# lines(density(BETA.store[,2],adj=2),col="gray",lwd=2)
# legend( "topright",legend=c( expression(theta[mileage]),expression((beta)[mileage])), 
#         lwd=c(2,2),col=c("black","gray"),bty="n") 
# 
# plot(density(THETA.b[,7],adj=2),xlim=range(BETA.store[,7]), 
#      main="",xlab="beta on mileage",ylab="posterior density",lwd=2)
# lines(density(BETA.store[,7],adj=2),col="gray",lwd=2)
# legend( "topright",legend=c( expression(theta[age]),expression((beta)[age])), 
#         lwd=c(2,2),col=c("black","gray"),bty="n") 

# Density of posteriors
par(mfrow = c(3,3), mar = c(5,4,3,1))
plot(density(BETA.ps[,1],adj=2),xlim=range(BETA.store[,1]),
     main="",xlab="Posterior mean of beta Intersection",ylab="Density",lwd=2)
plot(density(BETA.ps[,2],adj=2),xlim=range(BETA.store[,2]),
     main="",xlab="Posterior mean of beta Mileage",ylab="Density",lwd=2)
plot(density(BETA.ps[,3],adj=2),xlim=range(BETA.store[,3]),
     main="",xlab="Posterior mean of beta E85",ylab="Density",lwd=2)
plot(density(BETA.ps[,4],adj=2),xlim=range(BETA.store[,4]),
     main="",xlab="Posterior mean of beta Gasoline",ylab="Density",lwd=2)
plot(density(BETA.ps[,5],adj=2),xlim=range(BETA.store[,5]),
     main="",xlab="Posterior mean of beta Hybrid",ylab="Density",lwd=2)
plot(density(BETA.ps[,6],adj=2),xlim=range(BETA.store[,6]),
     main="",xlab="Posterior mean of beta Plug-In Hybrid",ylab="Density",lwd=2)
plot(density(BETA.ps[,7],adj=2),xlim=range(BETA.store[,7]),
     main="",xlab="Posterior mean of beta Age",ylab="Density",lwd=2)
plot(density(BETA.ps[,8],adj=2),xlim=range(BETA.store[,8]),
     main="",xlab="Posterior mean of beta Accident",ylab="Density",lwd=2)


# Visualization on shrinkage:
new_BETA.LS <- NULL
for (j in 1:m){
  fit <- lm(Y[[j]]~-1+X[[j]])
  new_BETA.LS <- rbind(new_BETA.LS, c(fit$coef))
}

par(mfrow = c(3,3), mar = c(5,4,3,1))

coef_names <- colnames(BETA.ps)

for(k in 1:(ncol(BETA.ps))){
  x <- N
  plot(N, new_BETA.LS[, k], pch = 16, col = "blue",
       xlab = "Group size", ylab = "Coefficient Value",
       main = coef_names[k], xaxt = "n")
  
  axis(1, at = x, las = 2, cex.axis = 0.6)
  points(N, BETA.ps[, k], pch = 17, col = "red")
  abline(h = beta_hat[k], col = "red", lwd = 2)
}
plot.new()
legend(
  "center",
  legend = c("OLS", "Bayesian"),
  col = c("blue", "red"),
  pch = c(16, 17),
  bty = "n"
)

# Statistics of posteriors
# -----------------------------
# CIs for posteriors:
sigma2_ci <- quantile(S2.b, probs = c(0.025, 0.975))
sigma2_ci
theta1_ci <- quantile(THETA.b[,1], probs = c(0.025, 0.975))
theta1_ci
theta2_ci <- quantile(THETA.b[,2], probs = c(0.025, 0.975))
theta2_ci
theta3_ci <- quantile(THETA.b[,3], probs = c(0.025, 0.975))
theta3_ci
theta4_ci <- quantile(THETA.b[,4], probs = c(0.025, 0.975))
theta4_ci
theta5_ci <- quantile(THETA.b[,5], probs = c(0.025, 0.975))
theta5_ci
theta6_ci <- quantile(THETA.b[,6], probs = c(0.025, 0.975))
theta6_ci
theta7_ci <- quantile(THETA.b[,7], probs = c(0.025, 0.975))
theta7_ci
theta8_ci <- quantile(THETA.b[,8], probs = c(0.025, 0.975))
theta8_ci
# means for posteriors
sigma2_mean <- mean(S2.b)
sigma2_mean
theta1_mean <- mean(THETA.b[,1])
theta1_mean
theta2_mean <- mean(THETA.b[,2])
theta2_mean
theta3_mean <- mean(THETA.b[,3])
theta3_mean
theta4_mean <- mean(THETA.b[,4])
theta4_mean
theta5_mean <- mean(THETA.b[,5])
theta5_mean
theta6_mean <- mean(THETA.b[,6])
theta6_mean
theta7_mean <- mean(THETA.b[,7])
theta7_mean
theta8_mean <- mean(THETA.b[,8])
theta8_mean
