# B_simulation_study.R
# Task 2 Part B:

# Set working directory to the script's location for reproducibility.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (!is.null(sys.function(1))) {
    ofile <- sys.frame(1)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  return(getwd())
}
setwd(get_script_dir())
# This is a simulation for the two-component Gaussian mixture model.
# This file makes two synthetic datasets:
# One has two bell curves far apart (clearly separated).
# One has overlapping bell curves (substantially overlapping).
# Runs the Gibbs sampler several times on each dataset, starting from different places.
# Produces trace plots, density plots, autocorrelation plots, and a few numeric summaries.
sink("output/simulation/console_output.txt")
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

# Outputs:
# coda needs an mcmc.list.
to_mcmc_list <- function(fits) {
  coda::mcmc.list(lapply(fits, function(f) coda::mcmc(f$samples)))
}
mcmc_sep <- to_mcmc_list(fit_sep)
mcmc_over <- to_mcmc_list(fit_over)

# Trace plots.
png("output/simulation/figures/trace_separated.png", width = 900, height = 600)
plot(mcmc_sep)
dev.off()
png("output/simulation/figures/trace_overlapping.png", width = 900, height = 600)
plot(mcmc_over)
dev.off()
# Posterior densities estimated from the draws.
d1 <- density(fit_sep[[1]]$samples[, "mu1"])
d2 <- density(fit_sep[[1]]$samples[, "mu2"])
png("output/simulation/figures/density_separated.png", width = 900, height = 600)
plot(d1, col = "blue", lwd = 2,
     main = "Posterior of mu1 and mu2 (separated)", xlab = "value",
     xlim = c(-4, 4), ylim = c(0, max(d1$y, d2$y)))
lines(d2, col = "red", lwd = 2)
legend("topright", legend = c("mu1", "mu2"), col = c("blue", "red"), lwd = 2)
dev.off()
d1 <- density(fit_over[[1]]$samples[, "mu1"])
d2 <- density(fit_over[[1]]$samples[, "mu2"])
png("output/simulation/figures/density_overlapping.png", width = 900, height = 600)
plot(d1, col = "blue", lwd = 2,
     main = "Posterior of mu1 and mu2 (overlapping)", xlab = "value",
     ylim = c(0, max(d1$y, d2$y)))
lines(d2, col = "red", lwd = 2)
legend("topright", legend = c("mu1", "mu2"), col = c("blue", "red"), lwd = 2)
dev.off()
# Autocorrelation measures how much a draw depends on the previous draw.
# High autocorrelation means the new draw is strongly influenced by the previous draw and doesn't give that much new information (overlapping components).
# With a flat likelihood, the new draw is heavily influenced by the previous one.
# In the separated case, since the data is clear, the sampler can jump to the correct answer independently without needing to rely on the previous draw.
png("output/simulation/figures/acf_separated.png", width = 900, height = 600)
autocorr.plot(mcmc_sep, auto.layout = TRUE)
dev.off()
png("output/simulation/figures/acf_overlapping.png", width = 900, height = 600)
autocorr.plot(mcmc_over, auto.layout = TRUE)
dev.off()
# Gelman-Rubin diagnostic compares variation within each chain to variation between chains.
# Values near 1 suggest the chains agree.
# Values above 1 suggest they are not gathered in the same region.
cat("\nGelman-Rubin (separated)\n")
print(gelman.diag(mcmc_sep))
cat("\nGelman-Rubin (overlapping)\n")
print(gelman.diag(mcmc_over))
# Effective sample size is how many independent samples the chain is worth.
# Lower ESS means noisier, less reliable posterior summaries (overlapping scenario).
cat("\nEffective sample size (separated)\n")
print(effectiveSize(mcmc_sep))
cat("\nEffective sample size (overlapping)\n")
print(effectiveSize(mcmc_over))
# Posterior summaries compute means and 95% credible intervals for each parameter.
# Labels can switch, so mu1 and mu2 may swap meaning across chains.
# Part C discusses this in detail, here the raw summaries are reported.
summarise <- function(fits, label) {
  draws <- do.call(rbind, lapply(fits, function(f) f$samples))
  means <- colMeans(draws)
  ci <- apply(draws, 2, quantile, probs = c(0.025, 0.975))
  cat("\nPosterior summary (", label, ")\n", sep = "")
  print(round(rbind(mean = means, ci), 3))
}
summarise(fit_sep, "separated")
summarise(fit_over, "overlapping")

# Component inference overlap issue:
# For the separated and overlapping cases, mu1 is arranged to be smaller than mu2.
# Then the posterior probability is computed.
# Reorder the pooled Z so component 1 is always the lower-mean component.
pool_ordered_Z <- function(fits) {
  z_all <- do.call(rbind, lapply(fits, function(f) f$Z)) # Z rows are iterations, columns are observations, values are 1 or 2.
  d_all <- do.call(rbind, lapply(fits, function(f) f$samples)) # d rows are iterations, columns are mu1, mu2, pi, sigma2.
  swap <- d_all[, "mu1"] > d_all[, "mu2"] # rows with mu1 > mu2 need to be swapped
  z_all[swap, ] <- 3L - z_all[swap, ] # 3 - 1 = 2, 3 - 2 = 1 flips labels in the Z table so component 1 is left mean, and component 2 is right mean.
  z_all
}
Z_sep <- pool_ordered_Z(fit_sep)
Z_over <- pool_ordered_Z(fit_over)
p1_sep <- colMeans(Z_sep == 1L) # Proportion of iterations in which each observation was assigned to component 1.
p1_over <- colMeans(Z_over == 1L)
# Fraction of observations with posterior probability between 0.1 and 0.9.
uncertain_sep <- mean(p1_sep > 0.1 & p1_sep < 0.9)
uncertain_over <- mean(p1_over > 0.1 & p1_over < 0.9)
cat("\nAllocation uncertainty (fraction with P(component 1) in (0.1, 0.9)):\n")
cat("separated:", round(uncertain_sep, 3), "\n")
cat("overlapping:", round(uncertain_over, 3), "\n")
# Histograms of allocation probabilities.
png("output/simulation/figures/allocation_separated.png", width = 900, height = 600)
hist(p1_sep, breaks = 20, col = "lightblue", border = "white",
     main = "Separated allocation probabilities",
     xlab = "posterior P(component 1)")
dev.off()
png("output/simulation/figures/allocation_overlapping.png", width = 900, height = 600)
hist(p1_over, breaks = 20, col = "lightpink", border = "white",
     main = "Overlapping allocation probabilities",
     xlab = "posterior P(component 1)")
dev.off()
cat("\nRun complete, figures are in output/simulation/figures/.\n")
sink()
