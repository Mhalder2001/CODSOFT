#1.#time#
library(datasets) 
data=Nile
data 
plot(data)
require(graphics)
data1<-ts(data,start=c(1871,1), end=c(1970,12), frequency=12) 
data1
decomp<-decompose(data1) # to decompose the data into different components plot(decomp)
Autocorr<-acf(data1,lag.max = 5,plot=TRUE) 
Autocorr
Partial_Autocorr<-acf(data1,lag.max = 5,type=c("partial"),plot=TRUE) 
Partial_Autocorr
#-----------------------------------------------
#2.
library(stats)
ts_data <- ts(c(13, 8, 15, 4, 4, 12, 11, 7, 14, 12))
autocorrelation <- function(data, lag = 3) 
{
acf_result <- acf(data, lag.max = lag, plot = FALSE)$acf 
return(acf_result)
}
partial_autocorrelation <- function(data, lag = 3)
{ 
pacf_result <- pacf(data, lag.max = lag, plot = FALSE)$acf 
return(pacf_result)
}
mean_ts <- mean(ts_data)
variance_ts <- var(ts_data)
autocorr_ts <- autocorrelation(ts_data)
partial_autocorr_ts <- partial_autocorrelation(ts_data)
cat("Mean of the time series data:", mean_ts, "\n")
cat("Variance of the time series data:", variance_ts, "\n") 
cat("Autocorrelation up to lag 3:", autocorr_ts, "\n")
cat("Partial autocorrelation up to lag 3:", partial_autocorr_ts, "\n")
#--------------------------------------------
#3.#multivariate#
head_length_first_son <- c(191, 195, 181, 176, 208, 189, 188, 192, 179, 183, 190, 188, 163, 186, 181, 192)
head_breadth_first_son <- c(155, 149, 148, 144, 157, 150, 152, 152, 158, 147, 159, 151, 137, 153, 140, 154)
head_length_second_son <- c(179, 201, 185, 171, 192, 190, 197, 186, 187, 174, 195, 187, 161, 173, 182, 185)
head_breadth_second_son <- c(145, 152, 149, 142, 152, 149, 159, 151, 148, 147, 157, 158, 130, 148, 146, 152)
data <- cbind(
  x1 = head_length_first_son,
  x2 = head_breadth_first_son,
  x3 = head_length_second_son,
  x4 = head_breadth_second_son
)
data
# 1. Mean Vector (MLE of μ)
mean_vector <-rowMeans(data) 
print(mean_vector)

# 2. Covariance Matrix (MLE of Σ)
cov_matrix <-cov(data)
cat("Covariance Matrix (MLE of Σ):\n")
print(cov_matrix)

# 3. Correlation Matrix (ρ)
cor_matrix <- cor(data)
print(cor_matrix)

# 4. Conditional Distribution Parameters
S11 <- cov_matrix[3:4, 3:4]      # Sub-matrix for (x1, x2)
S12 <- cov_matrix[3:4, 1:2]      # Sub-matrix between (x1, x2) and (x3, x4)
S21 <- t(S12)                    # Transpose of S12
S22 <- cov_matrix[1:2, 1:2]      # Sub-matrix for (x3, x4)
# Solving for S11 Inverse using Manual Inversion (2x2 matrix)
inv_S11<-solve(S11)
# Conditional Covariance: 
S22_1<-S21%*%inv_S11
S22.1 <- S22 - S21 %*% inv_S11 %*% S12
cat("Conditional Covariance Matrix (S22.1):\n")
print(S22.1)

# 5. Partial Correlation Between x3 and x4 Given x1, x2
r34<-cor_matrix[3,4]
r31<-cor_matrix[3,1]
r41<-cor_matrix[4,1]
r32<-cor_matrix[3,2]
r42<-cor_matrix[4,2]
numrator<-r34-(r31*r41)-(r32*r42)
denomenator<-sqrt((1-r31^2)*(1-r41^2)*(1-r32^2)*(1-r42^2))
r34.12<-numrator/denomenator
cat("Partial Correlation r34.12:\n", r34.12, "\n")

# 6. Fisher's Z for Confidence Interval
z_value <- 0.5 * log((1 + r34.12) / (1 - r34.12)) # Fisher's Z-transform
se <- 1 / sqrt(n - 3)
lower <- tanh(z_value - 1.96 * se)
upper <- tanh(z_value + 1.96 * se)
cat("95% Confidence Interval for r34.12:\n")
cat("Lower Bound:", lower, "\nUpper Bound:", upper, "\n")

# 7. Sample Multiple Correlation Coefficient
# Manual Calculation for x3 ~ (x1, x2)
model_x3<-lm(head_length_second_son~head_length_first_son +head_breadth_first_son)
R3_12<-summary(model_x3)$r.squared
R3_12
cat("Multiple Correlation Coefficient (x3 ~ x1, x2):\n", R3_12, "\n")
#Calculation for x4 ~ (x1, x2)
model_x4<-lm(head_breadth_second_son~head_length_first_son +head_breadth_first_son)
R4_12<-summary(model_x4)$r.squared
R4_12
cat("Multiple Correlation Coefficient (x3 ~ x1, x2):\n", R4_12, "\n")
#------------------------------machine learning--------
#1.kde
data <- c(5.65746599, 5.38283914, 2.79892121, 2.85423660, 2.95252721, 5.42626667,
          7.66239113, -0.18001073, 0.65083500, 2.40276530, -0.09929884, 6.32619215,
          5.03650752, 2.07470777, 1.78019174, 6.12891558, 4.05352439, 2.02686971,
          3.50834853, -2.76449768, 4.98428763, 3.01292677, 2.82448038, 3.98110437,
          5.09371862, 5.97961648, 4.56968496, -0.48814532, 5.08736697, 2.41757609)

# Kernel Density Estimate
kde <- density(data)
kde

# Naïve Density Estimator
naive_density <- function(x, data, h) {
  n <- length(data)
  return(sum((abs(x - data) <= h) / (2 * h)) / n)
}

# Naïve Density Estimator for a range of x
h <- 1.0 # Bandwidth for Naïve Density Estimation
x_range <- seq(min(data) - 1, max(data) + 1, length.out = 100)
naive_estimates <- sapply(x_range, naive_density, data = data, h = h)

# Plotting both
plot(kde, main = "Kernel Density Estimate vs Naïve Density Estimator", 
     col = "blue", lwd = 2, ylim = c(0, max(kde$y, naive_estimates)), xlab = "X", ylab = "Density")
lines(x_range, naive_estimates, col = "red", lwd = 2)
legend("topright", legend = c("KDE", "Naïve Density"), col = c("blue", "red"), lty = 1, lwd = 2)
####-----------------------------
#logistics
library(MASS)
data('birthwt')
birthwt
colnames(birthwt)
total_rows<-nrow(birthwt)
train_index<-sample(seq_len(total_rows),size=0.7*total_rows)
train_set<-birthwt[train_index,]
test_set<-birthwt[-train_index,]
train_set
#fit logistics 
model<-glm(low~age+lwt+race+smoke+ptl+ht+ui+ftv+bwt,data= train_set)
summary(model)
prediction<-predict(model,newdata=test_set,type="response")
predict_val<-ifelse(prediction>0.5,1,0)
predict_val
confussion_matrix<-table(Predicted=predict_val,Actual=test_set$low)
print(confussion_matrix)
TN<-confussion_matrix[1,1]
FP<-confussion_matrix[1,2]
FN<-confussion_matrix[2,1]
TP<-confussion_matrix[2,2]
sensitivity<-TP/(TP+FN)
specificity<-TN/(TN+FP)
ppv<-TP/(TP+FP)
npv<-TN/(TN+FN)
sensitivity
specificity
ppv
npv
#--------------------------------------------BIOSTAT----
# Data for city A (relaxed gun laws)
incidents_A <- 10
population_A <- 100000
# Data for city B (strict gun laws)
incidents_B <- 5
population_B <- 100000
# Calculate incidence rates
incidence_A <- incidents_A / population_A
incidence_B <- incidents_B / population_B
# Calculate Relative Risk (RR)
relative_risk_A <- incidence_A / incidence_B
relative_risk_A
relative_risk_B <- incidence_B / incidence_A
relative_risk_B
# Calculate Risk Difference (RD)
risk_difference <- incidence_A - incidence_B
# Output the results for part 1

cat("Risk Difference (RD) =", risk_difference, "\n\n")
#c)The seemingly obvious conclusion is that the relaxed gun laws in city. A cause more gun violence, quintupling the risk. However, before jumping to conclusions, it may be helpful to consider the following questions:
•Is the age distribution and socioeconomic status of each Population similar? Younger people involved in gangs or individuals of low socio economic status, may more likely to resort to gun violence. City A may be more prone to such situations.
•Were the risk exposure patterns several decades ago, when the laws were first induced, similar to those in the present? "Are the judicial systems and records of gun violenes, different in each city?

#-----------------------------------------------------------
# PART 2: Breast Cancer and Calcium Supplements
# Data for cases and controls
cases_no_calcium <- 75
cases_total <- 100
controls_no_calcium <- 25
controls_total <- 100

# Create a table to display the data
data_table <- matrix(c(cases_no_calcium, cases_total - cases_no_calcium, 
                       controls_no_calcium, controls_total - controls_no_calcium), 
                     nrow = 2, byrow = TRUE)
colnames(data_table) <- c("No Calcium", "Calcium")
rownames(data_table) <- c("Cases", "Controls")

# Display the data table
cat("Data Table:\n")
print(data_table)

# Calculate odds for cases and controls
odds_cases <- cases_no_calcium / (cases_total - cases_no_calcium)
odds_controls <- controls_no_calcium / (controls_total - controls_no_calcium)

# Calculate Odds Ratio (OR)
odds_ratio <- odds_cases / odds_controls

# Calculate prevalence of non-calcium intake in cases and controls
prevalence_cases <- cases_no_calcium / cases_total
prevalence_controls <- controls_no_calcium / controls_total

# Output the results for part 2
cat("\nOdds of exposure in cases =", odds_cases, "\n")
cat("Odds of exposure in controls =", odds_controls, "\n")
cat("Odds Ratio (OR) =", odds_ratio, "\n")
cat("Prevalence in cases =", prevalence_cases, "\n")
cat("Prevalence in controls =", prevalence_controls)
###d>After calculating the odds ratio, we observe a 3-Fold differences in the prevalence rate (75% vs 25%) change to a 9 - Fold differences in the odds ratio. Clearly, the two methods produce opposing results.
#----------------------------------------------------
# Given data are:
non_obese_non_exposed <- 0.03
non_obese_exposed <- non_obese_non_exposed * 1.2
obese_non_exposed <- non_obese_non_exposed * (2 / 3)
obese_exposed <- obese_non_exposed * 1.2

# Create a data frame to display incidence rates
incidence_table <- data.frame(
  Group = c("Non-Obese, Non-Exposed", "Non-Obese, Exposed", "Obese, Non-Exposed", "Obese, Exposed"),
  Incidence_Rate = c(non_obese_non_exposed, non_obese_exposed, obese_non_exposed, obese_exposed)
)
print(incidence_table)

# Calculating Odds and Odds Ratio for Non-Obese group as an example
odds_non_exposed <- non_obese_non_exposed / (1 - non_obese_non_exposed)
odds_exposed <- non_obese_exposed / (1 - non_obese_exposed)
odds_ratio <- odds_exposed / odds_non_exposed

# Displaying Odds Ratio and Relative Risk
cat("Odds Ratio (Non-Obese):", odds_ratio, "\n")
cat("Relative Risk (given): 1.2\n")
###(c)Overall, we can see that decreasing the baseline incidence will decrease the odd ratio (3.00 in those who are non-obese versus 1.23 in these who are obese). Obviously, these
#------------------------------------------
---------------------------------------------------
no_folate_neural <- 631            # Neural tube defects (No folate)
folate_neural <- 24                     # Neural tube defects (Folate)
no_folate_premature <- 727    # Premature births (No folate)
folate_premature <- 563          # Premature births (Folate)

# Attributable Risk (AR) calculations
AR_neural <- no_folate_neural - folate_neural
AR_premature <- no_folate_premature - folate_premature

cat("Attributable Risk for Neural Tube Defects:", AR_neural, "per 100,000\n")
cat("Attributable Risk for Premature Births:", AR_premature, "per 100,000\n")
##So, we claim of pregnant women not consuming folate, 96.2% of neural tube defect cases can be attributed to a lack of folate supplementation. Therefore, if the cause were to be removed, the disease could be reduced by up to 96.2% and 607 lives could be saved. Similarly, the attributable risk for premature births is 22.6%.
------------------------------------------------------------------
###-----------------------LDA--
# Given data
lifetimes <- c(22.3, 26.8, 30.3, 31.9, 32.1, 33.3, 33.7, 33.9, 34.7, 36.1, 36.4, 36.5, 
               36.6, 37.1, 37.6, 38.2, 38.5, 38.7, 38.7, 38.9, 38.9, 39.1, 41.1, 41.1, 
               41.4, 42.4, 43.6, 43.8, 44.0, 45.3, 45.8, 50.4, 51.3, 51.4, 51.5)
# Log-transform the data
log_lifetimes <- log(lifetimes)

# Estimate parameters (mu and sigma) using MLE
mu_hat <- mean(log_lifetimes)
sigma_hat <- sd(log_lifetimes)

# Calculate Mean and Median Time to Failure
mean_time_to_failure <- exp(mu_hat + (sigma_hat^2 / 2))
median_time_to_failure <- exp(mu_hat)

# Survival Function
survival_function <- function(t) {
  1 - pnorm((log(t) - mu_hat) / sigma_hat)
}
# Hazard Function
hazard_function <- function(t) {
  dnorm((log(t) - mu_hat) / sigma_hat) / (t * survival_function(t))
}
# Generate a sequence of times for plotting
time_seq <- seq(min(lifetimes), max(lifetimes), length.out = 100)
# Survival and Hazard values
survival_vals <- survival_function(time_seq)
hazard_vals <- hazard_function(time_seq)

# Plot Survival Curve
plot(time_seq, survival_vals, type = "l", col = "blue", lwd = 2,
     xlab = "Time (weeks)", ylab = "Survival Probability",
     main = "Survival Curve")
# Plot Hazard Curve
plot(time_seq, hazard_vals, type = "l", col = "red", lwd = 2,
     xlab = "Time (weeks)", ylab = "Hazard Rate",
     main = "Hazard Curve")
# results
cat("Estimated Parameters:\n")
cat("Mu (mean of log-lifetimes):", mu_hat, "\n")
cat("Sigma (std. dev. of log-lifetimes):", sigma_hat, "\n\n")
cat("Mean Time to Failure (MTTF):", mean_time_to_failure, "weeks\n")
cat("Median Time to Failure:", median_time_to_failure, "weeks\n")
-------------------------------------------------------
#100 SAMPLE WEIBULL 
shape_param <- 3     # Shape 
scale_param <- 10    # Scale
n_samples <- 100     # Number of observations

# Weibull-distributed data
data <- scale_param * (-log(runif(n_samples)))^(1 / shape_param)

# Negative log-likelihood function
neg_log_likelihood <- function(alpha, beta) {
  if (alpha <= 0 || beta <= 0) return(Inf) # Ensure parameters are positive
  n <- length(data)
  log_likelihood <- 
    n * log(alpha) - n * log(beta) +
    (alpha - 1) * sum(log(data)) -
    sum((data / beta)^alpha)
  return(-log_likelihood)  # Negative log-likelihood
}

# Generate 2D grid for shape (alpha) and scale (beta)
alpha_vals <- seq(2, 4, length.out = 50)  # Shape parameter grid
beta_vals <- seq(8, 12, length.out = 50)  # Scale parameter grid
likelihood_matrix <- matrix(0, nrow = length(alpha_vals), ncol = length(beta_vals))

# Computing log-likelihood for each combination
for (i in 1:length(alpha_vals)) {
  for (j in 1:length(beta_vals)) {
    likelihood_matrix[i, j] <- -neg_log_likelihood(alpha_vals[i], beta_vals[j])
  }}
# 2D likelihood plot
filled.contour(
  x = alpha_vals, y = beta_vals, z = likelihood_matrix,
  xlab = "Shape (alpha)", ylab = "Scale (beta)",
  main = "2D Likelihood Plot of Weibull Model"
)
# MLE Estimation
result <- optim(par = c(2, 5), fn = function(params) neg_log_likelihood(params[1], params[2]), 
                method = "L-BFGS-B", lower = c(1e-5, 1e-5))
mle_shape <- result$par[1]
mle_scale <- result$par[2]

# Mean failure time
mean_failure_time <- mle_scale * gamma(1 + 1 / mle_shape)  # Using Weibull formula
sample_mean <- mean(data)  # Sample mean
#results
cat("MLE Shape (alpha):", mle_shape, "\n")
cat("MLE Scale (beta):", mle_scale, "\n")
cat("MLE Mean Failure Time:", mean_failure_time, "\n")
cat("Sample Mean Failure Time:", sample_mean, "\n")
------------------------------------------------------------------------------
#LEUKAEMIA
# Given data
lifetimes <- c(22.3, 26.8, 30.3, 31.9, 32.1, 33.3, 33.7, 33.9, 34.7, 36.1, 36.4, 36.5, 
               36.6, 37.1, 37.6, 38.2, 38.5, 38.7, 38.7, 38.9, 38.9, 39.1, 41.1, 41.1, 
               41.4, 42.4, 43.6, 43.8, 44.0, 45.3, 45.8, 50.4, 51.3, 51.4, 51.5)
# Log-transform the data
log_lifetimes <- log(lifetimes)

# Estimate parameters (mu and sigma) using MLE
mu_hat <- mean(log_lifetimes)
sigma_hat <- sd(log_lifetimes)

# Calculate Mean and Median Time to Failure
mean_time_to_failure <- exp(mu_hat + (sigma_hat^2 / 2))
median_time_to_failure <- exp(mu_hat)

# Survival Function
survival_function <- function(t) {
  1 - pnorm((log(t) - mu_hat) / sigma_hat)
}
# Hazard Function
hazard_function <- function(t) {
  dnorm((log(t) - mu_hat) / sigma_hat) / (t * survival_function(t))
}
# Generate a sequence of times for plotting
time_seq <- seq(min(lifetimes), max(lifetimes), length.out = 100)
# Survival and Hazard values
survival_vals <- survival_function(time_seq)
hazard_vals <- hazard_function(time_seq)

# Plot Survival Curve
plot(time_seq, survival_vals, type = "l", col = "blue", lwd = 2,
     xlab = "Time (weeks)", ylab = "Survival Probability",
     main = "Survival Curve")
# Plot Hazard Curve
plot(time_seq, hazard_vals, type = "l", col = "red", lwd = 2,
     xlab = "Time (weeks)", ylab = "Hazard Rate",
     main = "Hazard Curve")
# results
cat("Estimated Parameters:\n")
cat("Mu (mean of log-lifetimes):", mu_hat, "\n")
cat("Sigma (std. dev. of log-lifetimes):", sigma_hat, "\n\n")
cat("Mean Time to Failure (MTTF):", mean_time_to_failure, "weeks\n")
cat("Median Time to Failure:", median_time_to_failure, "weeks\n")
---------------------extra-------------



