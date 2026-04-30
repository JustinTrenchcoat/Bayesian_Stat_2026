# this is the Bayesian page
# in here we would use grouped data to fit multiple OLS
# but we still use ungrouped data to create groups in Bayesian
# we will use info from multiple OLS to build our Bayesian Prior
# -------------------------------
# necessary library
# -------------------------------
library(coda)
# ------------------------------
# Load data
# -------------------------------
data <- read.csv("data/grouped.csv")
# -------------------------------
# Data preprocessing
# -------------------------------
data$fuel_type <- factor(data$fuel_type)
data$brand     <- factor(data$brand)
data$accident  <- factor(data$accident)
# -------------------------------
m_mean <- mean(data$mileage)
m_sd   <- sd(data$mileage)

a_mean <- mean(data$age)
a_sd   <- sd(data$age)

data$mileage <- (data$mileage - m_mean) / m_sd
data$age <- (data$age - a_mean) / a_sd

# -------------------------------
# Group setup
# -------------------------------
ids <- levels(data$brand)
m   <- length(ids)
Y <- vector("list", m)
X <- vector("list", m)
N <- numeric(m)

fit_list <- vector("list", m)
names(fit_list) <- ids

for(j in 1:m){
  
  idx <- which(data$brand == ids[j])
  
  group_dat <- droplevels(data[idx, ])
  N[j] <- nrow(group_dat)
  terms <- c("mileage", "age")
  
  if(length(unique(group_dat$fuel_type)) > 1){
    terms <- c(terms, "fuel_type")
  }
  
  if(length(unique(group_dat$accident)) > 1){
    terms <- c(terms, "accident")
  }
  
  form <- as.formula(
    paste("price ~", paste(terms, collapse = " + "))
  )
  
  fit_list[[j]] <- lm(form, data = group_dat)
}
fit <- lm(price ~ mileage + fuel_type + age + accident, data = data)
summary(fit)

# -------------------------------
# Extract coefficients
# -------------------------------
coef_names <- names(coef(fit))

coef_mat <- matrix(NA, nrow = m, ncol = length(coef_names))
colnames(coef_mat) <- coef_names
rownames(coef_mat) <- ids

for(j in 1:m){
  
  cf <- coef(fit_list[[j]])
  
  matched <- intersect(names(cf), coef_names)
  
  coef_mat[j, matched] <- cf[matched]
}

# -------------------------------
# View coefficient matrix
# -------------------------------
print(coef_mat)

# -------------------------------
# Plot beta vs group size
# -------------------------------
par(mfrow = c(3,3), mar = c(4,4,3,1))

for(k in colnames(coef_mat)){
  
  plot(
    N,
    coef_mat[,k],
    pch = 19,
    col = "blue",
    xlab = "Group Size",
    ylab = "Coefficient",
    main = k
  )
  
  abline(
    h = coef(fit)[k],
    col = "red",
    lwd = 2
  )
}

# -------------------------------
# Optional labeled plot for mileage_z
# -------------------------------
par(mfrow = c(1,1))

plot(
  N,
  coef_mat[, "mileage"],
  pch = 19,
  col = "blue",
  xlab = "Group Size",
  ylab = "Mileage Coefficient",
  main = "Mileage Coefficient vs Group Size"
)

text(
  N,
  coef_mat[, "mileage"],
  labels = ids,
  pos = 4,
  cex = 0.6
)

abline(
  h = coef(fit)["mileage"],
  col = "red",
  lwd = 2
)

#### Hierarchical regression model
# new dataset
ids<-sort(unique(data$brand)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-data[data$brand==ids[j], 6] 
  N[j]<- sum(data$brand==ids[j])
  x1<-data[data$brand==ids[j], 1]
  x2<-data[data$brand==ids[j], 2] 
  x3<-data[data$brand==ids[j], 3] 
  x4<-data[data$brand==ids[j], 4] 
  
  X[[j]]<-cbind( rep(1,N[j]), x1, x2, x3, x4)
}

## mvnormal simulation
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
# hyperparameter setup
# -----------------------------
BETA_clean <- BETA.LS[complete.cases(BETA.LS), ]
p<-dim(X[[1]])[2]

theta<-mu0<-apply(BETA_clean,2,mean)
BETA <- matrix(NA, m,p)

s2<-s20<-mean(S2.LS)
nu0<-1
eta0<-p+2 
Sigma<-S0<-L0<-cov(BETA_clean) 
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
BETA.ps  <- BETA*0
BETA.store <- list()

# -----------------------------
# MCMC
# -----------------------------
set.seed(632)

nsim   <- 50000
thin   <- 10

for(s in 1:nsim){
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
    Ej<-Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[[j]]/s2 )
    BETA[j,]<-rmvnorm(1,Ej,Vj) 
  } 
  ##
  
  ##update theta
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum))
  theta<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta <- matrix(theta, m, p, byrow = TRUE)
  SS <- S0 + t(BETA - mtheta) %*% (BETA - mtheta)
  iSigma <- rwish(1, eta0 + m, solve(SS))
  ##
  
  ##update s2
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

  ##store results
  if(s%%10==0) 
  { 
    cat(s,s2,"\n")
    S2.b<-c(S2.b,s2);THETA.b<-rbind(THETA.b,t(theta))
    Sigma.ps<-Sigma.ps+solve(iSigma) ; BETA.ps<-BETA.ps+BETA
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
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
effectiveSize(BETA.pp[,2])
effectiveSize(BETA.pp[,3])
effectiveSize(BETA.pp[,4])
effectiveSize(BETA.pp[,5])
effectiveSize(BETA.pp[,6])
effectiveSize(BETA.pp[,7])
effectiveSize(BETA.pp[,8])



# acf checks
# acf(S2.b)
# acf(THETA.b[,1])
# acf(BETA.pp[,1])

# Visualization of posterior coefficients
# customize brand is subject to change to any brand in the data set:
j <- which(ids=="Acura")
# customize coefficient search, two variables are subject to change:
k1 <- which(colnames(BETA.store[[1]]) == "mileage")
k2 <- which(colnames(BETA.store[[1]]) == "age")
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
# test shrinkage:
# ---------------------------------------
# Posterior group betas vs group size
# One plot per coefficient
# Blue points = posterior mean by group
# Red line   = global pooled OLS coefficient
# ---------------------------------------

par(mfrow = c(3,3), mar = c(4,4,3,1))

coef_names <- colnames(BETA.ps)

for(k in 1:ncol(BETA.ps)){
  plot(
    N,
    BETA.ps[, k],
    pch = 19,
    col = "blue",
    xlab = "Group Size",
    ylab = "Posterior Mean",
    main = coef_names[k]
  )
  
  # text(
  #   N,
  #   BETA.ps[, k],
  #   labels = ids,
  #   pos = 4,
  #   cex = 0.55
  # )
  
  abline(
    h = beta_hat[k],
    col = "red",
    lwd = 2
  )
}
