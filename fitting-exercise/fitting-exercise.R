## ---- packages --------
library(tidymodels)
library(readr)
library(dplyr)
library(ggplot2)

## ---- load-data --------
rawdat <- readr::read_csv("Mavoglurant_A2121_nmpk.csv")

## ---- explore-data --------
summary(rawdat)
skimr::skim(rawdat)
p1 <- rawdat %>% ggplot() +
  geom_line( aes( x = TIME, y = DV, group = as.factor(ID), color = as.factor(DOSE)) ) +
  facet_wrap( ~ DOSE, scales = "free_y")
plot(p1)

## ---- process-data --------
# remove those with OCC = 2
dat1 <- rawdat %>% filter(OCC == 1)

# total drug as a sum
dat_y <-  dat1 %>% filter((AMT == 0)) %>%
  group_by(ID) %>% 
  dplyr::summarize(Y = sum(DV)) 

#keep only time = 0 entry, it contains all we need
dat_t0 <- dat1 %>% filter(TIME == 0)

# merge data  
dat_merge <- left_join(dat_y, dat_t0, by = "ID")

# keep only useful variables
# also convert SEX and RACE to factors
dat <- dat_merge %>% select(Y,DOSE,AGE,SEX,RACE,WT,HT) %>% mutate(across(c(SEX, RACE), as.factor)) 
readr::write_rds(dat,"mavoglurant.rds")


# fit the linear models with Y as outcome 
# first model has only DOSE as predictor
# second model has all variables as predictors
lin_mod <- linear_reg() %>% set_engine("lm")
linfit1 <- lin_mod %>% fit(Y ~ DOSE, data = dat)
linfit2 <- lin_mod %>% fit(Y ~ ., data = dat)

# Compute the RMSE and R squared for model 1
metrics_1 <- linfit1 %>% 
  predict(dat) %>% 
  bind_cols(dat) %>% 
  metrics(truth = Y, estimate = .pred)

# Compute the RMSE and R squared for model 2
metrics_2 <- linfit2 %>% 
  predict(dat) %>% 
  bind_cols(dat) %>% 
  metrics(truth = Y, estimate = .pred)

# Print the results
print(metrics_1)
print(metrics_2)


## ---- fit-data-logistic --------
# fit the logistic models with SEX as outcome 
# first model has only DOSE as predictor
# second model has all variables as predictors
log_mod <- logistic_reg() %>% set_engine("glm")
logfit1 <- log_mod %>% fit(SEX ~ DOSE, data = dat)
logfit2 <- log_mod %>% fit(SEX ~ ., data = dat)

# Compute the accuracy and AUC for model 1
m1_acc <- logfit1 %>% 
  predict(dat) %>% 
  bind_cols(dat) %>% 
  metrics(truth = SEX, estimate = .pred_class) %>% 
  filter(.metric == "accuracy") 
m1_auc <-  logfit1 %>%
  predict(dat, type = "prob") %>%
  bind_cols(dat) %>%
  roc_auc(truth = SEX, .pred_1)


# Compute the accuracy and AUC for model 2
m2_acc <- logfit2 %>% 
  predict(dat) %>% 
  bind_cols(dat) %>% 
  metrics(truth = SEX, estimate = .pred_class) %>% 
  filter(.metric %in% c("accuracy"))
m2_auc <-  logfit2 %>%
  predict(dat, type = "prob") %>%
  bind_cols(dat) %>%
  roc_auc(truth = SEX, .pred_1)

# Print the results
print(m1_acc)
print(m2_acc)
print(m1_auc)
print(m2_auc)


## ---- plot-model-fit --------


# Augment linear regressions into data frame
mavo_lm3_aug <- augment(mavo_lm3, new_data = train_data) %>%
  select("Y", ".pred", ".resid") # select only the Y column and predictions
names(mavo_lm3_aug)[names(mavo_lm3_aug) == '.pred'] <- 'lm3_pred' # change name of prediction variable
names(mavo_lm3_aug)[names(mavo_lm3_aug) == '.resid'] <- 'lm3_resid' # change name of residual variable

mavo_lm4_aug <- augment(mavo_lm4, new_data = train_data)%>%
  select("Y", ".pred", ".resid")
names(mavo_lm4_aug)[names(mavo_lm4_aug) == '.pred'] <- 'lm4_pred'
names(mavo_lm4_aug)[names(mavo_lm4_aug) == '.resid'] <- 'lm4_resid'

mavo_lm0_aug <- augment(mavo_lm0, new_data = train_data) %>%
  select("Y", ".pred", ".resid") # select y and prediction column
names(mavo_lm0_aug)[names(mavo_lm0_aug) == '.pred'] <- 'lm0_pred'
names(mavo_lm0_aug)[names(mavo_lm0_aug) == '.resid'] <- 'lm0_resid'

combined_pred <- mavo_lm3_aug %>%
  left_join(mavo_lm4_aug, by='Y') %>%  
  left_join(mavo_lm0_aug, by='Y') # left join all three augmented data frames by Y
str(combined_pred) # Check the structure of the combined data frame

# Plot the combined data frame
ggplot(combined_pred, aes(x = Y)) + # add observed Y values
  geom_point(aes(y = lm3_pred, color = "Dose v. Y model")) + # add aes layer for lm3 predictions with a color for model name
  geom_point(aes(y = lm4_pred, color = "Full model")) + # add aes layer for lm4 predictions with a color for model name
  geom_point(aes(y = lm0_pred, color = "Null model")) + # add aes layer for lm0 predictions with a color for model name
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + # create 45 degree line
  labs(x = "Outcome (Y)", y = "Model Predictions", title = "Observed Value of Y versus Model Predictions") + # add axis labels and title
  scale_color_manual(values = c("Dose v. Y model" = "red", "Full model" = "blue", "Null model" = "green")) + # add colors for each model name
  theme_minimal()

## ---- residual-plot --------

# Residual plot
ggplot(mavo_lm4_aug, aes(x = lm4_pred, y = lm4_resid)) +
  geom_point() +  # Scatter plot of residuals
  geom_hline(yintercept = 0, linetype = "dashed") +  # add a dashed horizontal line at y = 0
  labs(x = "Predicted Values", y = "Residuals", title = "Residual Plot") +
  coord_cartesian(ylim = c(-2500, 2500)) +  # set y-axis range to -2500, 2500
  theme_minimal()

## ---- bootstrap --------

# set seed
set.seed(rngseed)

# bootstrapping
lm4_boot <- bootstraps(train_data, times = 100)  # bootstrap the train data 100 times

boot_pred <- array(NA, dim = c(nrow(train_data), length(lm4_boot))) # create an array vector to store the bootstrap data

# creating a loop
for (i in 1:length(lm4_boot)) {
  
  bootstrap <- analysis(lm4_boot$splits[[i]]) # extract the split samples from the bootstrap 
  model <- lm(Y ~ ., data = bootstrap) # fit the model for the bootstrap using full model formula and bootstrap sample data
  pred <- predict(model, newdata = train_data) # predict the new data from training data
  boot_pred[,i] <- pred # put the predictions into the list vector
}

dat_bs <- analysis(lm4_boot$splits[[i]]) # extract individual bootstrap samples
print(dat_bs)

# compute the mean and confidence intervals 
mean_ci <- function(x) { # use functions(x) to list the functions to calculate both the mean and CI
  mean <- mean(x) #mean(x) applies mean function to unspecified vector
  ci <- quantile(x, c(0.025, 0.5, 0.975)) #applies 95% CI quantiles to vector
  return(c(mean, ci))
}

preds <- apply(boot_pred, 2, mean_ci) # apply the customized function for mean and ci

preds_t <- t(preds) # t() transposes the matrix into a data frame

plot_data <- data.frame( #make a new dataframe to contain the observed Y, predictions from bootstrap, and bounds
  Y = c(train_data$Y),  # Observed values
  point_estimates = c(mavo_lm4_aug$lm4_pred),  # Point estimates from full model
  median = preds_t[, 2],  # Median from bootstrap sampling
  lower_CI = preds_t[, 3],  # Lower confidence limit from bootstrap sampling
  upper_CI = preds_t[, 4]   # Upper confidence limit from bootstrap sampling
)

# plot the point estimates and bounts along with the observed Y
ggplot(plot_data, aes(x = Y, y = point_estimates)) +
  geom_point(color = "black") +  # add layer for point estimates
  geom_point(aes(y = median), color = "blue") +  # add layer for median
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), color = "red", width = 0.2) +  # add layer for confidence intervals
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  # add 45 degree line
  labs(x = "Y Values Observed", y = "Mean Predicted Values", title = "Observed vs Predicted Values with Bootstrap Confidence Intervals") +  # Axes labels and title
  theme_minimal() 
