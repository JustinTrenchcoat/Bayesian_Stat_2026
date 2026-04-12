# Coding page!

traditional <- read.csv("data/traditional.csv")
green <- read.csv("data/green.csv")

# we will divide those data sets into 4 parts
# traditional part:
n <- nrow(traditional)
idx <- sample(n)

splits <- split(idx, rep(1:4, each = n/4))

Y_t1_test <- list()
X_t <- list()

for (j in 1:4) {
  Y_t1_test[[j]] <- traditional[splits[[j]], ]$price
  X_t[[j]] <- cbind(traditional[splits[[j]],]$milage, 
               traditional[splits[[j]],]$accident,
               traditional[splits[[j]],]$age)
}
test <- lm(Y_t1_test[[1]]~1+X_t[[1]])

Y_t1 <- traditional[splits[[1]], "price"]
X_t1 <- traditional[splits[[1]], c('milage', "accident", "age")]

Y_t2 <- traditional[splits[[2]], "price"]
X_t2 <- traditional[splits[[2]], c('milage', "accident", "age")]

Y_t3 <- traditional[splits[[3]], "price"]
X_t3 <- traditional[splits[[3]], c('milage', "accident", "age")]

Y_t4 <- traditional[splits[[4]], "price"]
X_t4 <- traditional[splits[[4]], c('milage', "accident", "age")]
#######################################
# green part:
n <- nrow(green)
idx <- sample(n)

splits <- split(idx, rep(1:4, each = n/4))
Y_g1 <- green[splits[[1]], "price"]
X_g1 <- model.matrix(price ~ ., data = green[splits[[1]], ])

Y_g2 <- green[splits[[2]], "price"]
X_g2 <- model.matrix(price ~ ., data = green[splits[[2]], ])

Y_g3 <- green[splits[[3]], "price"]
X_g3 <- model.matrix(price ~ ., data = green[splits[[3]], ])

Y_g4 <- green[splits[[4]], "price"]
X_g4 <- model.matrix(price ~ ., data = green[splits[[4]], ])

#### OLS fits(we need 8 of them!!)
S2_1.LS <- BETA_1.LS <- NULL

fit <- lm(Y_t1~1+X_t1)
BETA_1.LS <- rbind(BETA_1.LS, c(fit$coef))
S2_1.LS <- c(S2_1.LS, summary(fit)$sigma2)
fit <- lm.fit(x=X_g1, y = Y_g1)
BETA_1.LS <- rbind(BETA_1.LS, c(fit$coef))
S2_1.LS <- c(S2_1.LS, summary(fit)$sigma2)
