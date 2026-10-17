# A_gibbs_sampler.R
# Task 2 Part A:
# This file defines the function gibbs_mixture(), which will be called from other R files.

# In the couse assignment, a 2 component means (mu1 and mu2) Gaussian mixture model is given.
# mixing proportion = pi, variance = sigma2, hidden labels = Z_i (says which component each observation came from)

# SUMMARY OF THE STEPS:

# There are 2 components in the gaussian mixture model.
# We don't know what data points belong to which component.
# We use gibbs sampling algorithm method to estimate the parameters of these components.
# Afterward, we can calculate the probability that any specific data point belongs to component 1.

# In one run/chain of the gibbs sampler, these happen:
# Random guess or separation of the 2 components (just for the start.).
# Instead of starting random, the data could be split in half or k means clustering algorithm could be used to separate the data and put mu1 and mu2 to the centers of these 2 groups.
# Based on this initial guess, each data point in a single line are given a probability between 1 and 0, representing the probability of it belonging to component 1 instead of component 2.
# After this step, a random number is generated between 0 and 1 for each data point and compared with the calculated probability.
# If the random number is greater than the calculated probability of component 1, that data point is assigned to component 2.
# The calculated values of mu1, mu2 (centers of component 1 and 2), variance, ratio of the components against eachother are changed by getting a random number from normal, inverse gamma and beta  distribution.
# mu1, mu2, sigma (variance), pi (ratio of components)) are registered into a matrix, then the algorithm updates these values with the random generated ones explained above, and continues the same steps like this first iteration.

# Since the early iterations are influenced by the starting values, they are discarded.

gibbs_mixture <- function(y, # take numeric vector y
n_iter = 10000, # iteratioins to run
burn = 2000, # how many initial iterations to throw away
m0 = 0, s0 = 10, a = 1, b = 1, alpha0 = 1, beta0 = 1, # numbers that define the prior distributions (our beliefs before seeing data)
init = NULL, # supply own starting labels. if null, function chooses itself
store_Z = FALSE) # control if labels are saved at every iteration (if so, it uses more memory)
{ n <- length(y) # store sample size

# Choosing the starting labels:
  if (is.null(init)) { 
    km <- kmeans(y, centers = 2) # split the data into 2 clusters, assign each point a label 1 or 2
    Z <- km$cluster
  } else {
    Z <- init
  }
# Purely random start often gives one component with almost no observations, making the first few iterations unstable.
# K-means gives a reasonable starting partition.

# Starting values for the parameters (computed from the  starting labels):
  mu1 <- mean(y[Z == 1]) # mean of points labeled 1
  mu2 <- mean(y[Z == 2]) # mean of points labeled 2
  pi <- mean(Z == 1) # fraction labeled 1
  sigma2 <- var(y) # overall sample variance

# Storage
  samples <- matrix(NA_real_, nrow = n_iter, ncol = 4) # empty matrix, one row per iteration, 4 columns. fills as it goes
  colnames(samples) <- c("mu1", "mu2", "pi", "sigma2") # faster, cleaner than growing a vector each iteration

  if (store_Z) { # optional storage for the labels at every iteration
    Z_store <- matrix(NA_integer_, nrow = n_iter, ncol = n) 
  } else {
    Z_store <- NULL
  }
# Only will store the matrix in Part D where we need posterior uncertainty about allocations (store_Z = TRUE).
# Otherwise the n_iter by n matrix is big and uses a lot of memory. Default is "(store_Z = FALSE)"

