# Coding page!

data <- read.csv("data/cleaned.csv")

ids<-sort(unique(data$fuel_type)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-data[data$fuel_type==ids[j], "price"] 
  N[j]<- sum(data$fuel_type==ids[j])
  x1<-data[data$fuel_type==ids[j],"milage"] 
  x2<-data[data$fuel_type==ids[j],"accident"]
  x3<-data[data$fuel_type==ids[j],"age"]
  X[[j]]<-cbind( rep(1,N[j]), x1,x2,x3  )
}

#### OLS fits
S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[[j]]~-1+X[[j]] )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 
