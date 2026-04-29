library(rsample)
set.seed(632)

data <- read.csv("data/grouped.csv")
data$fuel_type <- factor(data$fuel_type)
data$brand <- factor(data$brand)
data$accident <- factor(data$accident)

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

pairs(data[,-5], pch=19, col='orange', lower.panel=panel.smooth)


# preliminary OLS fit
fit<-lm(price~1+mileage+fuel_type+age+accident, data = data_train)
summary(fit)
pred <- predict(fit, newdata = data_test)
plot(data_test$price, pred,
     xlab = "Actual Price",
     ylab = "Predicted Price",
     main = "Predicted vs Actual on OLS")
abline(0,1,col="red")
