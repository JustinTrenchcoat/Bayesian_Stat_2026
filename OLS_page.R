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

