### Code to compute ENT (effective number of tests) ###
#To obtain exposures correlation-adjusted p-value threshold

# 1. Compute the correlation matrix
# remove ID column, if any, and keep only tested exposures columns
cor_matrix <- cor(exposuresDF, use = "pairwise.complete.obs", method = "spearman")

# 2. Perform eigen decomposition
eigen_values <- eigen(cor_matrix)$values

# 3. Calculate ENT (effective number of tests) using the adapted formula
ENT <- sum((eigen_values - 1) * (eigen_values > 1))

# 4. Adjust the Bonferroni threshold
alpha <- 0.05  # Nominal significance level
adjusted_threshold <- alpha / ENT

# Output results
cat("Effective Number of Tests (ENT):", round(ENT, 3), "\n") #3.561
cat("Adjusted p-value threshold:", round(adjusted_threshold, 3), "\n") # 0.014
