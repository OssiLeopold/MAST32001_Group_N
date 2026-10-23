# C_label_switching.R
# Task 2 Part C: 

# The file looks at label switching and how to identify the two component Gaussian mixture.
# The likelihood and posterior are the same even if the labels are swithced.
# (mu1, mu2, pi)  <->  (mu2, mu1, 1 - pi)
# The posterior has at least two symmetric modes.
# A Gibbs sampler can move between them.
# When it does, mu1 and mu2 swap meaning across the chains.
# This is bad because if the labels are changed in the different runs, then the results will mix with each other.

# This script first shows the problem.
# Then fixes it with an after processing relabelling step in the start that makes mu1 < mu2.

# Setup:
sink("output/label_switching/console_output.txt")
source("A_gibbs_sampler.R")
set.seed(2026)

# Making fake data and run chains from opposite starts:
# 1000 numbers, 500 from a bell-shaped cloud around -2.5, 500 from a bell-shaped cloud around +2.5.
# The groups are moderately separated. They overlap a little, but not completely.
# If the groups are too far apart, labels almost never switch.
# If they overlap too much, it becomes hard to tell what is going on.
# A moderate separation shows the label-switching effect clearly.
# The components are identifiable, but the posterior still has two symmetric modes corresponding to the two labellings.

# A single chain often stays in one mode.
# To expose the two modes we run several chains from deliberately opposite starts and combine them.

n1 <- 500
y_mid <- c(rnorm(n1, mean = -2.5, sd = 1), rnorm(n1, mean = 2.5, sd = 1))

# Initial labelling putting component 1 on a chosen side.
init_side <- function(y, side = c("left", "right")) {
  side <- match.arg(side)
  if (side == "left") as.integer(y > 0) + 1L # 1 = left, 2 = right
  else as.integer(y < 0) + 1L # 1 = right, 2 = left
}

fits <- list(
  left1 = gibbs_mixture(y_mid, n_iter = 4000, burn = 1000, init = init_side(y_mid, "left"), store_Z = FALSE),
  left2 = gibbs_mixture(y_mid, n_iter = 4000, burn = 1000, init = init_side(y_mid, "left"), store_Z = FALSE),
  right1 = gibbs_mixture(y_mid, n_iter = 4000, burn = 1000, init = init_side(y_mid, "right"), store_Z = FALSE),
  right2 = gibbs_mixture(y_mid, n_iter = 4000, burn = 1000, init = init_side(y_mid, "right"), store_Z = FALSE)
)

# Combining all chains so they mix both labelling versions.
draws <- do.call(rbind, lapply(fits, function(f) f$samples))

# Trace plots of individual chains.
# If chain named left, component 1 is on the left, if chain right, component 1 is on the right.
# The chains are stable, but in 2 different chains, the labels are switched.
# A single chain would show only one mode and hide the label switching.
png("output/label_switching/figures/label_switching_trace.png", width = 900, height = 600)
par(mfrow = c(2, 2))
plot(fits$left1$samples[, "mu1"],  type = "l", col = "blue", main = "mu1, chain started left",  ylab = "value", xlab = "iteration")
plot(fits$right1$samples[, "mu1"], type = "l", col = "blue", main = "mu1, chain started right", ylab = "value", xlab = "iteration")
plot(fits$left1$samples[, "mu2"],  type = "l", col = "red", main = "mu2, chain started left",  ylab = "value", xlab = "iteration")
plot(fits$right1$samples[, "mu2"], type = "l", col = "red", main = "mu2, chain started right", ylab = "value", xlab = "iteration")
par(mfrow = c(1, 1))
dev.off()
# Joint scatter of (mu1, mu2).
# If there wasn't label switching, there would only be 1 clster.
png("output/label_switching/figures/label_switching_joint.png", width = 700, height = 700)
plot(draws[, "mu1"], draws[, "mu2"], pch = 20, cex = 0.3,
     col = rgb(0, 0, 0, 0.05),
     xlab = "mu1", ylab = "mu2",
     main = "Joint posterior of (mu1, mu2) with label switching")
abline(h = 0, lty = 2, col = "gray")
abline(v = 0, lty = 2, col = "gray")
dev.off()

# misleading summaries
# The averages of mu1 will be averaged with both the labelling groups, making it misleading.
summary_table <- function(m) {
  means <- colMeans(m)
  ci <- apply(m, 2, quantile, probs = c(0.025, 0.975))
  rbind(mean = means, `2.5%` = ci[1, ], `97.5%` = ci[2, ])
}
cat("\nlabel switching ignored summaries\n")
print(round(summary_table(draws), 3))

# These do not depend on which component is called 1 or 2:
# sigma2 (variance)
# |mu1 - mu2| distance between two means.
# min(mu1, mu2) and max(mu1, mu2)the hill locations stay the same, just getting taller

cat("\nnot depending on label ordering\n")
invariants <- cbind(
  sigma2 = draws[, "sigma2"],
  dist_means = abs(draws[, "mu1"] - draws[, "mu2"]),
  lower_mean = pmin(draws[, "mu1"], draws[, "mu2"]),
  upper_mean = pmax(draws[, "mu1"], draws[, "mu2"])
)
print(round(summary_table(invariants), 3))

# After processing to get rid of label switching:
# In every iteration where mu1 > mu2, swap the meants and make 1 - pi instead of pi.

swap <- draws[, "mu1"] > draws[, "mu2"]
draws_ordered <- draws
draws_ordered[swap, "mu1"] <- draws[swap, "mu2"]
draws_ordered[swap, "mu2"] <- draws[swap, "mu1"]
draws_ordered[swap, "pi"]  <- 1 - draws[swap, "pi"]

cat("\nfixed after processing\n")
print(round(summary_table(draws_ordered), 3))


# Visual comparison before and after ordering
# Before ordering, m1 is splited into two hills, after, it has only one.
png("output/label_switching/figures/label_switching_before_after.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
plot(density(draws[, "mu1"]), col = "blue", lwd = 2,
     main = "Before ordering", xlab = "value",
     xlim = range(c(draws[, "mu1"], draws[, "mu2"])))
lines(density(draws[, "mu2"]), col = "red", lwd = 2)
legend("topright", legend = c("mu1", "mu2"),
       col = c("blue", "red"), lwd = 2)
plot(density(draws_ordered[, "mu1"]), col = "blue", lwd = 2,
     main = "After ordering (mu1 < mu2)", xlab = "value",
     xlim = range(c(draws_ordered[, "mu1"], draws_ordered[, "mu2"])))
lines(density(draws_ordered[, "mu2"]), col = "red", lwd = 2)
legend("topright", legend = c("mu1", "mu2"),
       col = c("blue", "red"), lwd = 2)
par(mfrow = c(1, 1))
dev.off()

cat("\nLabel switching part complete. Figures written to output/label_switching/figures/.\n")

sink()