data <- read.csv("data/cleaned.csv")
data$fuel_type <- factor(data$fuel_type)
data$brand <- factor(data$brand)
data$accident <- factor(data$accident)

data_size<- nrow(data)
set.seed(632)
train_ind <- sample(seq_len(data_size), size = floor(0.8 * data_size))

data_train <- data[train_ind, ]
data_test <- data[-train_ind, ]

#standardizing the numerical variables
data_train$mileage <- scale(data_train$mileage)
data_train$age <- scale(data_train$age)

data_test$mileage <- scale(data_test$mileage)
data_test$age <- scale(data_test$age)
pairs(data_train, pch=19, col='orange', lower.panel=panel.smooth)


# preliminary OLS fit
fit<-lm(price~1+mileage+fuel_type+age+accident, data = data_train)
summary(fit)
# plot(as.numeric(fitted(fit)), as.numeric(residuals(fit)), type="p")

##############################################################################
# beyond this is purely experimental

# now if we fit them in groups:
ids<-sort(unique(data$brand)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
Y_test <- list()
X_test <- list()
N_test <- NULL
for(j in 1:m)
{
  Y[[j]]<-data_train[data_train$brand==ids[j], "price"]
  N[j]<- sum(data_train$brand==ids[j])
  x1<-data_train[data_train$brand==ids[j],"mileage"] 
  x2<-data_train[data_train$brand==ids[j],"fuel_type"]
  x3<-data_train[data_train$brand==ids[j],"age"]
  x4 <- data_train[data_train$brand==ids[j], "accident"]
  X[[j]]<-cbind( rep(1,N[j]), x1,x2,x3,x4)
  #################################
  Y_test[[j]]<-data_test[data_test$brand==ids[j], "price"]
  N_test[j]<- sum(data_test$brand==ids[j])
  x1_test<-data_test[data_test$brand==ids[j],"mileage"] 
  x2_test<-data_test[data_test$brand==ids[j],"fuel_type"]
  x3_test<-data_test[data_test$brand==ids[j],"age"]
  x4_test <- data_test[data_test$brand==ids[j], "accident"]
  X_test[[j]]<-cbind( rep(1,N_test[j]), x1_test,x2_test,x3_test, x4_test)
}
colnames(X[[j]]) <- c("Intercept","mileage","fuel_type","age", "accident")
colnames(X_test[[j]]) <- c("Intercept","mileage","fuel_type","age", "accident")


#### OLS fits
S2.LS<-BETA.LS<-NULL
OLS_fit_list <- list()
for(j in 1:m) {
  fit_group<-lm(Y[[j]]~ -1+X[[j]] )
  OLS_fit_list[[j]] <-fit_group
  BETA.LS<-rbind(BETA.LS,c(fit_group$coef)) 
  S2.LS<-c(S2.LS, summary(fit_group)$sigma^2) 
} 

model_ok <- sapply(OLS_fit_list, function(mod) {
  coef_valid <- !any(is.na(coef(mod)))
  sigma <- summary(mod)$sigma
  sigma_valid <- !is.na(sigma) && is.finite(sigma)
  
  coef_valid && sigma_valid
})

all(model_ok)

# tests:
RMSE_test <- NULL
RMSE_train <- NULL

for(j in 1:m){
  
  fit_test <- OLS_fit_list[[j]]
  
  # train predictions
  pred_train <- predict(fit_test)
  
  rmse_train <- sqrt(mean((Y[[j]] - pred_train)^2))
  
  # test predictions
  pred_test <- as.vector(X_test[[j]] %*% coef(fit_test))
  
  rmse_test <- sqrt(mean((Y_test[[j]] - pred_test)^2))
  
  RMSE_train <- c(RMSE_train, rmse_train)
  RMSE_test  <- c(RMSE_test, rmse_test)
}

data.frame(Group=1:m, Train_RMSE=RMSE_train, Test_RMSE=RMSE_test)

# age
plot( range(data_train$age),range(data_train$price),type="n",xlab="Car age", 
      ylab="log(Price)")
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,4],col="gray")  }
abline(fit$coefficients[1], fit$coefficients[4], col="black")
# mileage:
plot( range(data_train$mileage),range(data_train$price),type="n",xlab="Mileage", 
      ylab="log(Price)")
for(j in 1:m) {   abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }
abline(fit$coefficients[1], fit$coefficients[2], col="black")



BETA.MLS<-apply(BETA.LS,2,mean)
abline(BETA.MLS[1],BETA.MLS[4],lwd=2)

plot(N,BETA.LS[,1],xlab="sample size",ylab="intercept")
abline(h= BETA.MLS[1],col="black",lwd=2)
plot(N,BETA.LS[,4],xlab="sample size",ylab="slope_age")
abline(h= BETA.MLS[4],col="black",lwd=2)
