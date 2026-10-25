# D_Application.R
# Task 2 Part D: Faithful dataset

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
# This script fits a two-component Gaussian mixture to the Old Faithful eruption data.
# It relabels the posterior draws to keep mu1 < mu2, checks convergence,
# summarizes the parameters, and examines uncertainty in the component allocations.

# Setup:
source("A_gibbs_sampler.R")
library(coda)
set.seed(2026)

sink("output/faithful/console_output.txt")

# This loads the eruption-duration data into a variable
data(faithful)
y <- faithful$eruptions

# Run three chains from different initial partitions.
# The first one uses a standard 2-part kmeans split.
# The second one uses the same standard kmeans split with the labels reversed,
# in order to determine if results depend on initial labeling.
# Last one splits the observation by the median
km_cluster <- kmeans(y, centers = 2)$cluster

chains <- list(
  kmeans = gibbs_mixture(
    y, n_iter = 5000, burn = 1000,
    init = km_cluster, store_Z = TRUE
  ),
  kmeans_swapped = gibbs_mixture(
    y, n_iter = 5000, burn = 1000,
    init = 3L - km_cluster, store_Z = TRUE
  ),
  median = gibbs_mixture(
    y, n_iter = 5000, burn = 1000,
    init = as.integer(y > median(y)) + 1L, store_Z = TRUE
  )
)

# Checks whether the means of component 1 are larger than component 2
# If component 1 means are greater than component 2, the means, labels and related data is also swapped
for (k in seq_along(chains)) {
  samples <- chains[[k]]$samples
  swap <- samples[, "mu1"] > samples[, "mu2"]

  original_means <- samples[swap, c("mu1", "mu2"), drop = FALSE]
  samples[swap, "mu1"] <- original_means[, "mu2"]
  samples[swap, "mu2"] <- original_means[, "mu1"]
  samples[swap, "pi"] <- 1 - samples[swap, "pi"]

  chains[[k]]$samples <- samples
  chains[[k]]$Z[swap, ] <- 3L - chains[[k]]$Z[swap, ]
}

# Check convergence and effective sample sizes across chains.
# Takes the parameters that were sampled from each of the chains and formats it
# in a way that is usable by coda (coda mcmc object).
mcmc_runs <- mcmc.list(lapply(chains, function(fit) mcmc(fit$samples)))
cat("\nGelman-Rubin convergence diagnostics\n")
print(gelman.diag(mcmc_runs))
cat("\nEffective sample sizes\n")
print(effectiveSize(mcmc_runs))

# Combine the parameter draws and summarize them.
# For each parameter it calculated the posterior mean, the 2.5th percentile,
# and the 97.5th percentile.
all_draws <- do.call(rbind, lapply(chains, function(fit) fit$samples))
summary_table <- function(m) {
  means <- colMeans(m)
  ci <- apply(m, 2, quantile, probs = c(0.025, 0.975))
  rbind(mean = means, `2.5%` = ci[1, ], `97.5%` = ci[2, ])
}

# Prints all the summaries, rounded to three decimal points.
cat("\nPosterior summaries\n")
print(round(summary_table(all_draws), 3))

# Creates a two-panel png image showing fitted mixture and allocation uncertainty.
# Plots a histogram of the eruption durations, where it is scaled as a density in
# order to be compared with density curves.
# Means, weights and standard deviations are calculated and fitted
png("output/faithful/figures/faithful_mixture_and_allocations.png",
    width = 1100, height = 500)
par(mfrow = c(1, 2))

hist(y, breaks = 20, freq = FALSE, col = "lightgray",
     main = "Old Faithful fitted mixture", xlab = "Eruption duration (minutes)")
xg <- seq(min(y) - 0.5, max(y) + 0.5, length.out = 300)
m1 <- mean(all_draws[, "mu1"])
m2 <- mean(all_draws[, "mu2"])
w1 <- mean(all_draws[, "pi"])
sig <- sqrt(mean(all_draws[, "sigma2"]))

lines(xg, w1 * dnorm(xg, m1, sig) + (1 - w1) * dnorm(xg, m2, sig),
      col = "black", lwd = 2)
lines(xg, w1 * dnorm(xg, m1, sig), col = "blue", lty = 2)
lines(xg, (1 - w1) * dnorm(xg, m2, sig), col = "red", lty = 2)
legend("topright", legend = c("Mixture", "Component 1", "Component 2"),
       col = c("black", "blue", "red"), lty = c(1, 2, 2))

# Estimate each observation's posterior probability of belonging to component 1.
all_Z <- do.call(rbind, lapply(chains, function(fit) fit$Z))
p_comp1 <- colMeans(all_Z == 1L)

# Plot allocation probabilities against eruption duration; colors show the
# more likely component, and the dashed line marks the 0.5 cutoff.
plot(y, p_comp1, pch = 16,
     col = ifelse(p_comp1 > 0.5, "blue", "red"),
     xlab = "Eruption duration (minutes)", ylab = "P(Component 1)",
     main = "Allocation uncertainty")
abline(h = 0.5, lty = 2, col = "gray")

par(mfrow = c(1, 1))
dev.off()

cat("\nApplication complete. Graph saved in output/faithful/figures/.\n")
sink()