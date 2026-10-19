# B_simulation_study.R
# Task 2 Part B:

# This is a simulation for the two-component Gaussian mixture model.
# This file makes two synthetic datasets:
# One has two bell curves far apart (clearly separated).
# One has overlapping bell curves (substantially overlapping).
# Runs the Gibbs sampler several times on each dataset, starting from different places.
# Produces trace plots, density plots, autocorrelation plots, and a few numeric summaries.

source("A_gibbs_sampler.R")
library(coda)
set.seed(2026)

# Generate two datasets
# Situation 1: two bell curves far apart.
# Means at -3 and +3, same spread, equal mixing proportions.
# This is supposed to be easier as the two groups are almost visually distinguishable.
n1 <- 500
y_sep <- c(rnorm(n1, mean = -3, sd = 1), rnorm(n1, mean = 3, sd = 1))

# Situation 2 has two bell curves overlapping a lot.
# Means at -1 and +1, wider spread, equal mixing proportions.
# This is hard since the groups blend into one another.
y_over <- c(rnorm(n1, mean = -1, sd = 1.5), rnorm(n1, mean = 1, sd = 1.5))

# Run several Gibbs sampler chains on each dataset.
run_chains <- function(y, seeds) {
  fits <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    set.seed(seeds[i])
# Initial assignment of each data point with k-means clustering.
# Closer data put into same cluster, 2 groups in total. Better than random assignment.
# Different seeds give different k-means starts, so the chains do not all get stuck in the same wrong region.

# If the chain moves freely around the high posterior region, it is mixing well instead of getting stuck.
# When label switching occurs, both chains describe the same two groups but with mu1 and mu2 swapped.

# If the chains are exploring the same region, Gelman-Rubin is close to 1. (Situation 1, separated)
# If the chains are stuck in different regions, Gelman-Rubin is large for mu1 and mu2. (Situation 2, overlapping) (For Part C)
    fits[[i]] <- gibbs_mixture(
      y,
      n_iter = 6000,
      burn = 1000,
      store_Z = TRUE
    )
  }
  fits
}
seeds <- c(11, 22, 33, 44)
fit_sep <- run_chains(y_sep, seeds)
fit_over <- run_chains(y_over, seeds)
