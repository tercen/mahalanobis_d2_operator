suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr)
  library(tidyr)
})

ctx <- tercenCtx()

if (length(ctx$colors) < 1) stop("At least one color factor is required to define groups.")

# Get data as matrix: rows = observations (columns zone), cols = variables (rows zone)
mat <- ctx$as.matrix()
# mat has nrow = number of row factors (variables), ncol = number of column combinations (observations)
# Transpose so rows = observations, cols = variables
mat <- t(mat)

# Get group labels from color factors
# Build group vector from the .ci (column index) and color factors
n_obs <- nrow(mat)
n_var <- ncol(mat)

# Get the color factor values for each observation (column)
# Each column in the original crosstab is one observation
# We need to map .ci values to group labels
ci_values <- rep(0:(ncol(ctx$as.matrix()) - 1), each = 1)

# Get color factor data - colors are axis-level factors
color_df <- ctx$select(c(".ci", ctx$colors)) %>%
  distinct(.ci, .keep_all = TRUE) %>%
  arrange(.ci)

# Create composite group label
if (length(ctx$colors) == 1) {
  group_vec <- as.character(color_df[[ctx$colors[1]]])
} else {
  group_vec <- do.call(paste, c(color_df[ctx$colors], list(sep = ".")))
}

# Remove rows with NAs
complete_rows <- complete.cases(mat)
if (sum(complete_rows) < nrow(mat)) {
  mat <- mat[complete_rows, , drop = FALSE]
  group_vec <- group_vec[complete_rows]
}

groups <- unique(group_vec)
n_groups <- length(groups)
if (n_groups < 2) stop("At least 2 groups are required for Mahalanobis D².")

p <- ncol(mat)  # number of variables

# Compute group statistics
group_stats <- list()
for (g in groups) {
  idx <- which(group_vec == g)
  group_stats[[g]] <- list(
    n = length(idx),
    mean = colMeans(mat[idx, , drop = FALSE]),
    cov = cov(mat[idx, , drop = FALSE])
  )
}

# Compute pooled within-group covariance matrix
N <- nrow(mat)
S_pooled <- matrix(0, p, p)
for (g in groups) {
  n_g <- group_stats[[g]]$n
  S_pooled <- S_pooled + (n_g - 1) * group_stats[[g]]$cov
}
S_pooled <- S_pooled / (N - n_groups)

# Generate all unique pairs
pairs <- combn(groups, 2, simplify = FALSE)

# Compute D² for each pair and each method
results <- list()

for (pair in pairs) {
  g_i <- pair[1]
  g_j <- pair[2]

  mean_diff <- group_stats[[g_i]]$mean - group_stats[[g_j]]$mean

  # Method 1: Pooled covariance
  S_pooled_inv <- tryCatch(solve(S_pooled), error = function(e) {
    MASS::ginv(S_pooled)
  })
  d2_pooled <- as.numeric(t(mean_diff) %*% S_pooled_inv %*% mean_diff)

  # Method 2: Individual (average of two group covariances)
  S_avg <- (group_stats[[g_i]]$cov + group_stats[[g_j]]$cov) / 2
  S_avg_inv <- tryCatch(solve(S_avg), error = function(e) {
    MASS::ginv(S_avg)
  })
  d2_individual <- as.numeric(t(mean_diff) %*% S_avg_inv %*% mean_diff)

  # Method 3: Anderson (1973) asymptotic D² for unequal covariances
  # Anderson's formula: D² = delta' * [S_i/n_i + S_j/n_j]^{-1} * delta
  # This accounts for unequal covariance matrices and sample sizes
  n_i <- group_stats[[g_i]]$n
  n_j <- group_stats[[g_j]]$n
  S_anderson <- group_stats[[g_i]]$cov / n_i + group_stats[[g_j]]$cov / n_j
  S_anderson_inv <- tryCatch(solve(S_anderson), error = function(e) {
    MASS::ginv(S_anderson)
  })
  d2_anderson <- as.numeric(t(mean_diff) %*% S_anderson_inv %*% mean_diff)

  results <- c(results, list(
    tibble(group_i = g_i, group_j = g_j, d_squared = d2_pooled, method = "pooled"),
    tibble(group_i = g_i, group_j = g_j, d_squared = d2_individual, method = "individual"),
    tibble(group_i = g_i, group_j = g_j, d_squared = d2_anderson, method = "anderson")
  ))
}

result <- bind_rows(results)

# Add .ci and .ri for global summary output
result$.ci <- 0L
result$.ri <- 0L

result %>%
  ctx$addNamespace() %>%
  ctx$save()
