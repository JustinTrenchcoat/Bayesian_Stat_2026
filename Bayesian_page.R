# this is the Bayesian page
# in here we would use grouped data to fit multiple OLS
# but we still use ungrouped data to create groups in Bayesian
# we will use info from multiple OLS to build our Bayesian Prior
# -------------------------------
# necessary library
# -------------------------------
library(coda)
# -------------------------------
# load data
# -------------------------------
# data <- read.csv("data/cleaned.csv")
data_grouped <- read.csv("data/grouped.csv")

# -----------------------------
# Data preprocessing
# -----------------------------
data_grouped$fuel_type <- factor(data_grouped$fuel_type)
data_grouped$brand     <- factor(data_grouped$brand)
data_grouped$accident  <- factor(data_grouped$accident)

#standardizing the numerical variables
m_mean <- mean(data_grouped$mileage)
m_sd   <- sd(data_grouped$mileage)

a_mean <- mean(data_grouped$age)
a_sd   <- sd(data_grouped$age)


data_grouped$mileage <- (data_grouped$mileage - m_mean) / m_sd
data_grouped$age <- (data_grouped$age - a_mean) / a_sd

# -----------------------------
# Grouped OLS
# -----------------------------
ids<-sort(unique(data_grouped$brand)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-data_grouped[data_grouped$brand==ids[j], 6] 
  N[j]<- sum(data_grouped$brand==ids[j])
  x1<-data_grouped[data_grouped$brand==ids[j], 1]
  x2<-data_grouped[data_grouped$brand==ids[j], 2] 
  x3<-data_grouped[data_grouped$brand==ids[j], 3] 
  x4<-data_grouped[data_grouped$brand==ids[j], 4] 
  
  X[[j]]<-cbind( rep(1,N[j]), x1, x2, x3, x4)
}
#### OLS fits
S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[[j]]~-1+X[[j]] )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 

# OLS betas vs sample size
for (i in 1:4)

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
rmvnorm<-function(n,mu,Sigma)
{ 
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}

## Wishart simulation
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

## Setup
BETA_clean <- BETA.LS[complete.cases(BETA.LS), ]

p<-dim(X[[1]])[2]
theta<-mu0<-apply(BETA_clean,2,mean)
nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA_clean) ; BETA <- matrix(NA, m,p)
THETA.b<-S2.b<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps<-matrix(0,p,p)
SIGMA.PS<-NULL
BETA.ps<-BETA*0
BETA.pp<-NULL
set.seed(632)

## MCMC
for(s in 1:10000) {
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
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) ) )
  ##
  
  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]]%*%BETA[j,] )^2 ) }
  s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  ##
  ##store results
  if(s%%10==0) 
  { 
    cat(s,s2,"\n")
    S2.b<-c(S2.b,s2);THETA.b<-rbind(THETA.b,t(theta))
    Sigma.ps<-Sigma.ps+solve(iSigma) ; BETA.ps<-BETA.ps+BETA
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
  }
  ##
}



## MCMC diagnostics
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(THETA.b[,2])

apply(SIGMA.PS,2,effectiveSize)

tmp<-NULL;for(j in 1:dim(SIGMA.PS)[2]) { tmp<-c(tmp,acf(SIGMA.PS[,j])$acf[2]) }

acf(S2.b)
acf(THETA.b[,1])
acf(THETA.b[,2])
