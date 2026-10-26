# Load data
bodyfat <- read.csv('bodyfatmen.csv')



################## FULL MODEL WITH ALL REGRESSION VARIABLES ##################

bodyfat.model <- lm(density ~ ., data = bodyfat)

summary(bodyfat.model)

n <- nobs(bodyfat.model)              # sample size
p <- length(coef(bodyfat.model))      # number of parameters (incl. intercept)

# Density vs regressors
par(mfrow = c(4, 4))
plot(bodyfat$age, bodyfat$density)
plot(bodyfat$weight, bodyfat$density)
plot(bodyfat$height, bodyfat$density)
plot(bodyfat$neck, bodyfat$density)
plot(bodyfat$chest, bodyfat$density)
plot(bodyfat$abdomen, bodyfat$density)
plot(bodyfat$hip, bodyfat$density)
plot(bodyfat$thigh, bodyfat$density)
plot(bodyfat$knee, bodyfat$density)
plot(bodyfat$ankle, bodyfat$density)
plot(bodyfat$forearm, bodyfat$density)
plot(bodyfat$biceps, bodyfat$density)
plot(bodyfat$wrist, bodyfat$density)


######### Check Multicollinearity for full model #########
# Variable Inflation Factor for each regressor (cutoff value 5 or 10)
library(car)
vif(bodyfat.model)

library("GGally")  # ggpairs
ggpairs(data = bodyfat.model)


######### Detect Outliers #########

# Residuals vs. regressor variables
par(mfrow = c(4, 4))
plot(bodyfat$age, bodyfat.model$residuals)
plot(bodyfat$weight, bodyfat.model$residuals)
plot(bodyfat$height, bodyfat.model$residuals)
plot(bodyfat$neck, bodyfat.model$residuals)
plot(bodyfat$chest, bodyfat.model$residuals)
plot(bodyfat$abdomen, bodyfat.model$residuals)
plot(bodyfat$hip, bodyfat.model$residuals)
plot(bodyfat$thigh, bodyfat.model$residuals)
plot(bodyfat$knee, bodyfat.model$residuals)
plot(bodyfat$ankle, bodyfat.model$residuals)
plot(bodyfat$forearm, bodyfat.model$residuals)
plot(bodyfat$biceps, bodyfat.model$residuals)
plot(bodyfat$wrist, bodyfat.model$residuals)

# Cutoff values (from Belsley, Kuh, and Welsch)
leverage.cutoff <- 2 * p / n  # MPV p. 213
cooks.cutoff <- qf(0.50, p, n - p, lower.tail = FALSE)  # MPV p. 215
dfbetas.cutoff <- 2 / sqrt(n)  # MPV p. 218
dffits.cutoff <- 2 * sqrt(p / n)  # MPV p. 219
covratio_plus.cutoff <- 1 + 3*p/n
covratio_minus.cutoff <- 1 - 3*p/n
studres.cutoff <- qt(0.05 / 2, n - p, lower.tail = FALSE)  # MPV p. 135

# Residual plots
par(mfrow = c(1, 3))
plot(bodyfat.model, which = 1:3)

# Cook's distance plots
par(mfrow = c(1, 2))
plot(bodyfat.model, which = 4:5)

# DFBETAS and DFFITS, COVRATIO, COOK's DISTANCE
bodyfat.influence <- influence.measures(bodyfat.model)
summary(bodyfat.influence)


######### Evaluate Outliers #########

# CHECK HOW MUCH MS_RES DECREASES IN % WHEN REMOVING POINTS
mod39_bodyfat <- bodyfat[-c(39, 42), ]
mod39_bodyfat.model <- lm(density ~ ., data = mod39_bodyfat)
mod86_bodyfat <- bodyfat[-c(86, 42), ]
mod86_bodyfat.model <- lm(density ~ ., data = mod86_bodyfat)
mod96_bodyfat <- bodyfat[-c(96, 42), ]
mod96_bodyfat.model <- lm(density ~ ., data = mod96_bodyfat)

mod39_86_bodyfat <- bodyfat[-c(39, 42, 86), ]
mod39_86_bodyfat.model <- lm(density ~ ., data = mod39_86_bodyfat)

mod39_96_bodyfat <- bodyfat[-c(39, 42, 96), ]
mod39_96_bodyfat.model <- lm(density ~ ., data = mod39_96_bodyfat)

mod86_96_bodyfat <- bodyfat[-c(42, 86, 96), ]
mod86_96_bodyfat.model <- lm(density ~ ., data = mod86_96_bodyfat)

mod39_86_96_bodyfat <- bodyfat[-c(39, 42, 86, 96), ]
mod39_86_96_bodyfat.model <- lm(density ~ ., data = mod39_86_96_bodyfat)

residual_reduction_39 <- ((summary(bodyfat.model)$sigma^2 - summary(mod39_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100
residual_reduction_86 <- ((summary(bodyfat.model)$sigma^2 - summary(mod86_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100
residual_reduction_96 <- ((summary(bodyfat.model)$sigma^2 - summary(mod96_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100
residual_reduction_39_86 <- ((summary(bodyfat.model)$sigma^2 - summary(mod39_86_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100
residual_reduction_39_96 <- ((summary(bodyfat.model)$sigma^2 - summary(mod39_96_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100
residual_reduction_86_96 <- ((summary(bodyfat.model)$sigma^2 - summary(mod86_96_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100
residual_reduction_39_86_96 <- ((summary(bodyfat.model)$sigma^2 - summary(mod39_86_96_bodyfat.model)$sigma^2)/(summary(bodyfat.model)$sigma^2))*100

cat("Improvement removing 86:", round(residual_reduction_86, 2), "%\n")
cat("Improvement removing 39:", round(residual_reduction_39, 2), "%\n")
cat("Improvement removing 39 & 86:", round(residual_reduction_39_86, 2), "%\n")
cat("Improvement removing 96:", round(residual_reduction_96, 2), "%\n")
cat("Improvement removing 86 & 96:", round(residual_reduction_86_96, 2), "%\n")
cat("Improvement removing 39 & 96:", round(residual_reduction_39_96, 2), "%\n")
cat("Improvement removing 39 & 86 & 96:", round(residual_reduction_39_86_96, 2), "%\n")


######### Remove Outliers and see how it affects full model #########

# Residual plots with resplots.R (bootstrap) for model without outliers
source("resplot.R")
resplot(mod39_86_96_bodyfat.model)

# Partial residual plots
library(car)  # avPlots, boxcox, boxTidwell
avPlots(mod39_86_96_bodyfat.model)

# Studentized residuals
library(MASS)  # studres
par(mfrow = c(1, 2))
plot(studres(mod39_86_96_bodyfat.model))
plot(rstudent(mod39_86_96_bodyfat.model))

# Check Multicollinerity
library(car)
vif(mod39_86_96_bodyfat.model)

X <- model.matrix(mod39_86_96_bodyfat.model)
XtX <- cor(X[,-1])  # X^T X in correlation form

mod39_86_96_bodyfat.model.eigen <- eigen(XtX)
max(mod39_86_96_bodyfat.model.eigen$values) / min(mod39_86_96_bodyfat.model.eigen$values)

library("GGally")  # ggpairs
ggpairs(data = mod39_86_96_bodyfat.model)



################## SUBSET MODEL SELECTION ##################
# Stepwise regression (Variable selection)
step_mod <- step(mod39_86_96_bodyfat.model, direction = "both")
summary(step_mod)

library(leaps)

# Perform all possible regressions
all_models <- regsubsets(density ~ ., 
                         data = mod39_86_96_bodyfat,
                         nbest = 1,      # Keep best model of each size
                         nvmax = NULL)   # Consider all variables

# Get summary of regression
summary_models <- summary(all_models)

# Gather all relevant results (R-sqr, Adj R-sqr, Cp, BIC)
summary_models$which      # Which variables in each model
summary_models$rsq        # R-squared
summary_models$adjr2      # Adjusted R-squared
summary_models$cp         # Mallow's Cp
summary_models$bic        # Bayesian Information Criterion

# Create a comparison table
results <- data.frame(
  Num_Vars = 1:length(summary_models$rsq),
  R_squared = summary_models$rsq,
  Adj_R_squared = summary_models$adjr2,
  bic = summary_models$bic,
  Cp = summary_models$cp
)
print(results)

### Plot the submodelmodel results ###
par(mfrow = c(2, 2))

# Plot 1: R-squared
plot(results$Num_Vars, results$R_squared, 
     type = "b", pch = 19,
     xlab = "Number of Variables", 
     ylab = "R-squared",
     main = "R-squared vs Model Size")

# Plot 2: Adjusted R-squared
plot(results$Num_Vars, results$Adj_R_squared, 
     type = "b", pch = 19,
     xlab = "Number of Variables", 
     ylab = "Adjusted R-squared",
     main = "Adjusted R-squared vs Model Size")
points(which.max(results$Adj_R_squared), 
       max(results$Adj_R_squared), 
       col = "red", cex = 2, pch = 19)

# Plot 3: Mallow's Cp
plot(results$Num_Vars, results$Cp, 
     type = "b", pch = 19,
     xlab = "Number of Variables", 
     ylab = "Mallow's Cp",
     main = "Mallow's Cp vs Model Size")
abline(a = 1, b = 1, lty = 2, col = "blue")
legend("topright", legend = "Cp = p+1", lty = 2, col = "blue")

# Plot 4: BIC
plot(results$Num_Vars, results$bic, 
     type = "b", pch = 19,
     xlab = "Number of Variables", 
     ylab = "BIC",
     main = "BIC vs Model Size")
points(which.min(results$bic), 
       min(results$bic), 
       col = "green", cex = 2, pch = 19)

# See which variables are in each model 
coef(all_models, 3)
coef(all_models, 6)
coef(all_models, 7)
coef(all_models, 8)
coef(all_models, 9)

bodyfat6.model <- lm(density ~ age + weight + abdomen + thigh + forearm + wrist, 
                     data = mod39_86_96_bodyfat)
bodyfat7.model <- lm(density ~ age + weight + neck + abdomen + thigh + forearm + wrist, 
                     data = mod39_86_96_bodyfat)
bodyfat8.model <- lm(density ~ age + weight + neck + abdomen + hip + thigh + forearm + wrist, 
                     data = mod39_86_96_bodyfat)
###### PRESS ######
# PRESS calculation function
calc_press <- function(model) {
  pr <- residuals(model) / (1 - lm.influence(model)$hat)
  sum(pr^2)
}
library(MPV)

# Calculate PRESS for individual models
bodyfat6.press <- calc_press(bodyfat6.model)
bodyfat7.press <- calc_press(bodyfat7.model)
bodyfat8.press <- calc_press(bodyfat8.model)

# Baseline and full models
intercept.model <- lm(density ~ 1, data = mod39_86_96_bodyfat)  # Intercept only
baseline.model <- lm(density ~ ., data = mod39_86_96_bodyfat)  # All predictors

# Calculate PRESS for baseline models
intercept.press <- calc_press(intercept.model)
baseline.press <- calc_press(baseline.model)

# Compare all models
press.comparison <- data.frame(
  Model = c("Null (0 vars)", "Best 6-var", "Best 7-var", "Best 8-var", "All 13 vars"),
  PRESS = c(intercept.press, bodyfat6.press, bodyfat7.press, bodyfat8.press, baseline.press)
)

print(press.comparison)

# Select model 7 because of best results given above 



# Comparison of 7 regressor model vs full model (without outliers)
library(stargazer)

stargazer(mod39_86_96_bodyfat.model, bodyfat7.model,
          type = "text",
          title = "Model Comparison: Full Model vs Model 7",
          column.labels = c("Full Model", "Model 7"),
          keep.stat = c("n", "rsq", "adj.rsq", "ser", "f"))

######## Check multicollinearity for model with 7 regressors ########
# Look at VIFs for this model
vif(bodyfat7.model)
# VIF_weight > 10 --> weight is in multicollinearity with other regressors

# Check ggpairs matrix
library("GGally")  # ggpairs
ggpairs(data = bodyfat7.model)

# Check condition number
X <- model.matrix(bodyfat7.model)
XtX <- cor(X[,-1])  # X^T X in correlation form

bodyfat7.model.eigen <- eigen(XtX)
max(bodyfat7.model.eigen$values) / min(bodyfat7.model.eigen$values)
# No serious multicollinearity according to condition number

# We decide to remove "weight" because of its high VIF value


###### ASSEMBLE BEST SUBSET MODEL FROM REGSUBSETS WITHOUT WEIGHT ######

subset.model <- lm(density ~ age + neck + abdomen + thigh + forearm + wrist, 
              data = mod39_86_96_bodyfat)



################## ASSESS subset MODEL ##################

# ANOVA Summary of subset model
summary(subset.model)

# Residual plots for subset model
par(mfrow = c(1, 3))
plot(subset.model, which = 1:3)

# Cook's distance plots for subset model
par(mfrow = c(1, 2))
plot(subset.model, which = 4:5)

# VIFs for subset model
vif(subset.model) # No multicollinearity is seen

# Condition number for subset model
X <- model.matrix(subset.model)
XtX <- cor(X[,-1])  # X^T X in correlation form

subset.model.eigen <- eigen(XtX)
max(subset.model.eigen$values) / min(subset.model.eigen$values) # low condition number


#### Bootstrap for subset model ####
source("resplot.R")
resplot(subset.model)

library(car)
subset.case <- Boot(subset.model, method="case") #perform bootstrap resampling using the case method (data sampled with replacement)
subset.res <- Boot(subset.model, method="residual") #perform bootstrap resampling using the residual method (resample residuals)

confint(subset.case)
confint(subset.res)
summary(subset.model)
# We can see that thigh and forearm are insignificant

# Histogram of bootstrap estimates
par(mfrow = c(3, 2))
hist(subset.case)

par(mfrow = c(3, 2))
hist(subset.res)


######### ASSESS SMALLER MODELS WITHOUT THIGH AND FOREARM #########

smaller.model1 <- lm(density ~ age + neck + abdomen + wrist, 
                   data = mod39_86_96_bodyfat)
smaller.model2 <- lm(density ~ age + abdomen + wrist, 
                   data = mod39_86_96_bodyfat)

stargazer(subset.model, smaller.model1, smaller.model2,
          type = "text",
          title = "Model Comparison",
          column.labels = c("(age + neck + abdomen + thigh + forearm + wrist)", "(age + neck + abdomen + wrist)", "(age + abdomen + wrist)"),
          keep.stat = c("n", "rsq", "adj.rsq", "ser", "f"))
# We can see that for smaller.model1, neck is insignificant, for smaller.model2 all variables are significant


# Compare PRESS
original.press = calc_press(mod39_86_96_bodyfat.model)
subset.press = calc_press(subset.model)
smaller1.press = calc_press(smaller.model1)
smaller2.press = calc_press(smaller.model2)
print(original.press)
print(subset.press)
print(smaller1.press)
print(smaller2.press)

# Compare Mallow's Cp
library(olsrr)
ols_mallows_cp(mod39_86_96_bodyfat.model, mod39_86_96_bodyfat.model)
ols_mallows_cp(subset.model, mod39_86_96_bodyfat.model)
ols_mallows_cp(smaller.model1, mod39_86_96_bodyfat.model)
ols_mallows_cp(smaller.model2, mod39_86_96_bodyfat.model)
# Cp doesn't seem to change much between subset models

# Compare AIC and BIC
comparison.AIC <- AIC(mod39_86_96_bodyfat.model, subset.model, smaller.model1, smaller.model2)
print(comparison.AIC)

comparison.BIC <- BIC(mod39_86_96_bodyfat.model, subset.model, smaller.model1, smaller.model2)
print(comparison.BIC)

# Seems like the model with age + abdomen + wrist has about the same predictive
# power as larger models so we choose this model as our best by Occam's Razor. 



############# ASSESS FINAL MODEL (age + abdomen + wrist) #############

final.model <- smaller.model2

# Residual plots for final model
par(mfrow = c(2, 3))
plot(final.model, which = 1:5)

# Condition number for final model
X <- model.matrix(final.model)
XtX <- cor(X[,-1])  # X^T X in correlation form

final.model.eigen <- eigen(XtX)
max(final.model.eigen$values) / min(final.model.eigen$values) # very low condition number

# VIF for final model
vif(final.model) # very low VIFs

# Bootstrapping with final model
library(car)
final.case <- Boot(final.model, method="case") #perform bootstrap resampling using the case method (data sampled with replacement)
final.res <- Boot(final.model, method="residual") #perform bootstrap resampling using the residual method (resample residuals)

confint(final.case)
confint(final.res)
summary(final.model)

par(mfrow = c(3, 2))
hist(final.case)

par(mfrow = c(3, 2))
hist(final.res)

source("resplot.R")
resplot(final.model)
# No unusual deviation in residual plots


# Look at difference in performance for final model compared to original
final.residual_difference <- ((summary(mod39_86_96_bodyfat.model)$sigma^2 - summary(final.model)$sigma^2)/(summary(mod39_86_96_bodyfat.model)$sigma^2))*100
final.adjustedR2_difference <- ((summary(mod39_86_96_bodyfat.model)$adj.r.squared - summary(final.model)$adj.r.squared)/(summary(mod39_86_96_bodyfat.model)$adj.r.squared))*100
cat("Difference in residual between full model and final model:", round(final.residual_difference, 2), "%\n")
cat("Difference in adjusted R2 between full model and final model:", round(final.adjustedR2_difference, 2), "%\n")