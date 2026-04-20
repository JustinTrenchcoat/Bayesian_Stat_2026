data <- read.csv("data/cleaned.csv")
data_size<- nrow(data)
set.seed(632)
train_ind <- sample(seq_len(data_size), size = floor(0.8 * data_size))

data_train <- data[train_ind, ]
data_test <- data[-train_ind, ]

#standardizing the numerical variables
data_train$mileage <- scale(data_train$milage)
data_train$age <- scale(data_train$age)

data_test$mileage <- scale(data_test$milage)
data_test$age <- scale(data_test$age)

# now if we fit them in groups:
ids<-sort(unique(data$accident)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
Y_test <- list()
X_test <- list()
N_test <- NULL
for(j in 1:m)
{
  Y[[j]]<-data_train[data_train$accident==ids[j], "price"]
  N[j]<- sum(data_train$accident==ids[j])
  x1<-data_train[data_train$accident==ids[j],"mileage"] 
  x2<-data_train[data_train$accident==ids[j],"fuel_type"]
  x3<-data_train[data_train$accident==ids[j],"age"]
  X[[j]]<-cbind( rep(1,N[j]), x1,x2,x3)
  #################################
  Y_test[[j]]<-data_test[data_test$accident==ids[j], "price"]
  N_test[j]<- sum(data_test$accident==ids[j])
  x1_test<-data_test[data_test$accident==ids[j],"mileage"] 
  x2_test<-data_test[data_test$accident==ids[j],"fuel_type"]
  x3_test<-data_test[data_test$accident==ids[j],"age"]
  X_test[[j]]<-cbind( rep(1,N_test[j]), x1_test,x2_test,x3_test)
}
colnames(X[[j]]) <- c("Intercept","mileage","fuel_type","age")
colnames(X_test[[j]]) <- c("Intercept","mileage","fuel_type","age")

#### OLS fits
S2.LS<-BETA.LS<-NULL
OLS_fit_list <- list()
for(j in 1:m) {
  fit_group<-lm(Y[[j]]~ -1+X[[j]] )
  OLS_fit_list[[j]] <-fit_group
  BETA.LS<-rbind(BETA.LS,c(fit_group$coef)) 
  S2.LS<-c(S2.LS, summary(fit_group)$sigma^2) 
  print(summary(fit_group))
} 

#### Hierarchical regression model

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
p<-dim(X[[1]])[2]
theta<-mu0<-apply(BETA.LS,2,mean)
nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 
Sigma<-S0<-L0<-diag(diag(cov(BETA.LS))+1e-6)
BETA<-BETA.LS
THETA.b<-S2.b<-NULL
iL0<-solve(L0)  
iSigma<-solve(Sigma)
Sigma.ps<-matrix(0,p,p)
SIGMA.PS<-NULL
BETA.ps<-BETA*0
BETA.pp<-NULL
set.seed(1)
mu0[2]+c(-1.96,1.96)*sqrt(L0[2,2])

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
  BETA.store <- list()
  if(s%%10==0) 
  { 
    BETA.store[[length(BETA.store)+1]] <- BETA
    if (s%%1000==0){cat(s,s2,"\n")}
    S2.b<-c(S2.b,s2);THETA.b<-rbind(THETA.b,t(theta))
    Sigma.ps<-Sigma.ps+solve(iSigma) ; BETA.ps<-BETA.ps+BETA
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
  }
  ##
}


## MCMC diagnostics
library(coda)
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(THETA.b[,2])

apply(SIGMA.PS,2,effectiveSize)

tmp<-NULL;for(j in 1:dim(SIGMA.PS)[2]) { tmp<-c(tmp,acf(SIGMA.PS[,j])$acf[2]) }

acf(S2.b)
acf(THETA.b[,1])
acf(THETA.b[,2])

nkeep <- length(BETA.store)

RMSE_test <- NULL

for(j in 1:m){
  
  pred_mat <- matrix(0, nrow(X_test[[j]]), nkeep)
  
  for(s in 1:nkeep){
    beta_draw <- BETA.store[[s]][j, ]
    pred_mat[,s] <- X_test[[j]] %*% beta_draw
  }
  
  pred_mean <- rowMeans(pred_mat)
  
  rmse <- sqrt(mean((Y_test[[j]] - pred_mean)^2))
  
  RMSE_test <- c(RMSE_test, rmse)
}

RMSE_test


