# Author: Florian Schneider
# Date: June 3rd 2025

library(tidyverse)
library(ranger)
library(rsample)
library(foreach)
library(doMC)

# Concordance Correlation Coefficient (CCC)
ccc <- function(obs, pred) {
  mean_obs <- mean(obs)
  mean_pred <- mean(pred)
  var_obs <- var(obs)
  var_pred <- var(pred)
  cov_obs_pred <- cov(obs, pred)
  
  numerator <- 2 * cov_obs_pred
  denominator <- var_obs + var_pred + (mean_obs - mean_pred)^2
  
  return(numerator / denominator)
}

# R-squared
r2 <- function(obs, pred) {
  ss_total <- sum((obs - mean(obs))^2)
  ss_residual <- sum((pred - obs)^2)
  
  return(1 - (ss_residual / ss_total))
}

# RMSE
rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

# Load and preprocess data ####
dat <- read.csv(here(".../geodata.csv")) |>
  left_join(read.csv(here(".../EUROSTAT.csv"))) |>
  select(ID, all_of(readxl::read_xlsx(here(".../vars.xlsx")) |>
                      pull(name))) |>
  filter(OC < 10, RockFragment < 5)

# Load predictor variables
vars <- readxl::read_xlsx(here(".../vars.xlsx"))

# Filter data for grassland
grass <- dat |>
  filter(LU == "grassland") |>
  select(BD, all_of(vars$name[vars$grass == 1])) |>
  drop_na()

# Define tuning grid for ranger model
tuning_grid <- expand.grid(mtry = (1:27) * 2, 
                           min.node.size = (1:10) * 2)

# Number of outer and inner CV folds
k_outer <- 10 # Number of outer CV folds
k_inner <- 10 # Number of inner CV folds for tuning

# Nested-cross validation ####
# Define the outer CV folds
set.seed(123)
outer_folds <- vfold_cv(grass, v = k_outer)

# Outer loop for model evaluation
outer_test <- foreach(i = 1:k_outer, .combine = "rbind") %do% {
  # Status message
  print(paste(i, "of", k_outer, "outer folds"))
  
  # Split outer fold
  outer_train <- analysis(outer_folds$splits[[i]])
  outer_test <- assessment(outer_folds$splits[[i]])
  
  ## Inner loop ####
  # for trying out various models and their hyperparamters
  set.seed(123)
  inner_folds <- vfold_cv(outer_train, v = k_inner)
  
  cores <- parallel::detectCores() - 1
  registerDoMC(cores)
  
  inner_tuning <- foreach(j = 1:nrow(tuning_grid), .combine = "rbind") %:%
    foreach(k = 1:k_inner, .combine = "rbind") %dopar% {
      todo = tuning_grid[j, ]
      
      inner_train <- analysis(inner_folds$splits[[k]])
      inner_val <- assessment(inner_folds$splits[[k]])
      
      # Run model
      set.seed(123)
      m <- ranger(
        BD ~ .,
        data = inner_train,
        num.trees = 500,
        mtry = todo$mtry,
        min.node.size = todo$min.node.size
      )
      
      # Performance metrics from inner loop
      tibble(
        mtry = todo$mtry,
        min.node.size = todo$min.node.size,
        BD = inner_val$BD,
        BD_pred = as.numeric(predict(m, inner_val)$predictions)
      ) |>
        group_by(mtry, min.node.size) |>
        summarise(
          inner_ccc = ccc(obs = BD, pred = BD_pred),
          inner_r2 = r2(obs = BD, pred = BD_pred),
          inner_rmse = rmse(obs = BD, pred = BD_pred),
          inner_bias = mean(BD) - mean(BD_pred),
          inner_slope = coef(lm(BD_pred ~ BD))[2]
        ) |>
        ungroup()
      
    }
  
  try(registerDoSEQ(), silent = T)
  
  # Select best hyperparameters based on inner CV results
  inner_hyperparameters <- inner_tuning |>
    ungroup() |>
    arrange(-inner_ccc) |>
    slice(1)
  
  # Outer loop ####
  # To evaluate tuned model
  set.seed(123)
  m_tuned <- ranger(
    BD ~ .,
    data = outer_train,
    num.trees = 500,
    mtry = inner_hyperparameters$mtry,
    min.node.size = inner_hyperparameters$min.node.size
  )
  
  # Outer test result
  tibble(
    outer_fold = i,
    mtry = inner_hyperparameters$mtry,
    min.node.size = inner_hyperparameters$min.node.size,
    inner_ccc = inner_hyperparameters$inner_ccc,
    inner_slope = inner_hyperparameters$inner_slope,
    BD = as.numeric(outer_test$BD),
    BD_pred = predict(m_tuned, outer_test)$predictions
  )
  
}

# Summary of outer test results
outer_test_summary <- outer_test |>
  group_by(outer_fold, mtry, min.node.size, inner_ccc, inner_slope) |>
  summarise(
    outer_ccc = ccc(obs = BD, pred = BD_pred),
    outer_r2 = r2(obs = BD, pred = BD_pred),
    outer_rmse = rmse(obs = BD, pred = BD_pred),
    outer_bias = mean(BD) - mean(BD_pred),
    outer_slope = coef(lm(BD_pred ~ BD))[2]
  ) |>
  ungroup() |>
  pivot_longer(outer_ccc:outer_slope) |>
  group_by(name) |>
  summarise(mean = mean(value), sd = sd(value))

# Export results
writexl::write_xlsx(
  list(
    points = outer_test |>
      select(outer_fold, BD, BD_pred),
    metrics = outer_test_summary
  ),
  here(".../validation.xlsx")
)

# Production ####
cores = parallel::detectCores() - 1
registerDoMC(cores)

tuning <- foreach(j = 1:nrow(tuning_grid), .combine = "rbind") %:%
  foreach(k = 1:k_outer, .combine = "rbind") %dopar% {
    todo <- tuning_grid[j, ]
    
    train <- analysis(outer_folds$splits[[k]])
    val <- assessment(outer_folds$splits[[k]])
    
    # Run model
    set.seed(123)
    m <- ranger(
      BD ~ .,
      data = train,
      num.trees = 500,
      mtry = todo$mtry,
      min.node.size = todo$min.node.size
    )
    
    # Performance metrics
    tibble(
      mtry = todo$mtry,
      min.node.size = todo$min.node.size,
      BD = val$BD,
      BD_pred = as.numeric(predict(m, val)$predictions)
    ) |>
      group_by(mtry, min.node.size) |>
      summarise(
        ccc = ccc(obs = BD, pred = BD_pred),
        r2 = r2(obs = BD, pred = BD_pred),
        rmse = rmse(obs = BD, pred = BD_pred),
        bias = mean(BD) - mean(BD_pred),
        slope = coef(lm(BD_pred ~ BD))[2]
      ) |>
      ungroup()
    
  }

try(registerDoSEQ(), silent = T)

# Select best hyperparameters for production model
production_hyperparameters <- tuning |>
  ungroup() |>
  arrange(-ccc) |>
  slice(1)

set.seed(123)
m_production <- ranger(
  BD ~ .,
  data = grass,
  num.trees = 500,
  mtry = production_hyperparameters$mtry,
  min.node.size = production_hyperparameters$min.node.size,
  importance = "permutation"
)

# Export
saveRDS(m_production, here(".../grass-model.rds"))
