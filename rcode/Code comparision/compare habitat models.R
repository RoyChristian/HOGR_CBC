root <- getwd()

################################################################################
#                                 Packages
################################################################################
library(tidyverse)
library(reshape2)

################################################################################
#                                 Functions
################################################################################

# Function to load and join all cross-validation folds for a specific model
load_model_lppd <- function(model_name, folder_name, file_prefix, years = 41:47) {
  
  df_list <- lapply(years, function(y) {
    # Construct path based on directory structure
    file_path <- file.path(root, "output", folder_name, paste0(file_prefix, "_N", y, "lppd.csv"))
    
    # Read and keep only identifying strata and the lppd value
    read.csv(file_path) %>%
      select(Strata, lppd_pred) %>%
      rename(!!paste0("lppd_pred_", y) := lppd_pred)
  })
  
  # Join all years into a single wide dataframe 
  df_final <- df_list %>% 
    reduce(left_join, by = "Strata") %>%
    mutate(Model = model_name) %>%
    relocate(Strata, Model)
  
  return(df_final)
}

# Pairwise comparison function
lppd_zscore <- function(df_mod1, df_mod2) {
  
  # Extract numeric lppd columns (3 to end) and convert to single vectors
  # This aligns all 47 strata x 7 CV years into one vector of 329 values
  vec1 <- as.vector(as.matrix(df_mod1[, 3:ncol(df_mod1)]))
  vec2 <- as.vector(as.matrix(df_mod2[, 3:ncol(df_mod2)]))
  
  if(length(vec1) != length(vec2)) {
    stop("The two model dataframes must have the same dimensions for comparison.")
  }
  
  # Calculate pairwise differences (Model 1 - Model 2)
  diff <- vec1 - vec2
  
  # Calculate mean and standard error of the differences
  mean_diff <- mean(diff)
  se_diff <- sd(diff) / sqrt(length(diff))
  
  # Return the Z-score (Positive = Model 1 is better; Negative = Model 2 is better)
  return(mean_diff / se_diff)
}

################################################################################
#                       Load All Model Results
################################################################################

# Define model metadata for batch loading
model_meta <- list(
  list(id = "ponds_only",     fold = "ponds_only",     pre = "fit_breeding_ponds_only"),
  list(id = "noponds",        fold = "noponds",        pre = "fit_breeding_noponds"),
  list(id = "breeding_all",   fold = "breeding_all",   pre = "fit_breeding_all"),
  list(id = "winter",         fold = "winter",         pre = "fit_winter_only"),
  list(id = "noponds_winter", fold = "noponds_winter", pre = "fit_breeding_noponds_winter"),
  list(id = "linear_only",    fold = "linear_only",    pre = "fit_linear_only"),
  list(id = "climate_linear", fold = "climate_linear", pre = "fit_climate_linear")
)

# Apply the loading function to create a named list of dataframes
models <- lapply(model_meta, function(m) {
  load_model_lppd(m$id, m$fold, m$pre)
}) %>% setNames(sapply(model_meta, `[[`, "id"))

################################################################################
#                         Pairwise Comparisons
################################################################################

# We compare exactly two models at a time to determine relative predictive gain
comparisons <- list(
  # Habitat-specific comparisons
  list(a = "noponds",        b = "ponds_only"),
  list(a = "noponds",        b = "breeding_all"),
  
  # Breeding vs Winter drivers
  list(a = "breeding_all",   b = "winter"),
  list(a = "noponds",        b = "winter"),
  
  # Combined model performance
  list(a = "noponds_winter", b = "noponds"),
  list(a = "noponds_winter", b = "winter"),
  
  # Complexity vs Performance (Time/Linear components)
  list(a = "climate_linear", b = "noponds_winter"),
  list(a = "climate_linear", b = "linear_only")
)

# Run the Z-score calculation for each pair
results_list <- lapply(comparisons, function(pair) {
  z <- lppd_zscore(models[[pair$a]], models[[pair$b]])
  
  data.frame(
    Model_A = pair$a,
    Model_B = pair$b,
    Z_Score = round(z, 3),
    Best  = ifelse(z > 0, pair$a, pair$b)
  )
})

################################################################################
#                       Summary Table
################################################################################

final_comparison_table <- bind_rows(results_list)

print(final_comparison_table)

