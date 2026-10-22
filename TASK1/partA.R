library(JADE)

set.seed(20260816)

# Function to simulate autoregressive observations according to given ar and ma coefficients.
generate_arma = function(T, ar = numeric(0), ma = numeric(0)) {
  # Moving-average generation coefficients.
  psi = ARMAtoMA(ar = ar, ma = ma, lag.max = 2000)
  
  # Standard deviation for Gaussian innovation.
  sd_innov = 1 / sqrt(1 + sum(psi^2))
  
  # Generate the random autoregressive data.
  arima.sim(model = list(ar = ar, ma = ma), n = T, rand.gen = function(n) rnorm(n, sd = sd_innov))
}

# Function to generate signals for testing. T is the chosen sample size.
generate_signals = function(T) {
  s1 = generate_arma(T, ar = 0.6, ma = c(0, -0.4))
  s2 = generate_arma(T, ar = c(0.4, 0.1, 0.3), ma = c(0.2, 0.4))
  s3 = generate_arma(T, ar = 0.4)
  s4 = generate_arma(T, ma = c(0.5, 0.3, 0.1, 0.2))
  s5 = generate_arma(T, ar = c(0.1, 0.1, 0.3))
  
  signals = cbind(s1, s2, s3, s4, s5)
  
  return(signals)
}

# Function to generate noise according to setting 1: Gaussian white noise.
generate_noise = function(T) {
  noise = matrix(rnorm(T * 5), nrow = T, ncol = 5)
  
  return(noise)
}

# Function for generating random observations for testing.
generate_data = function(T) {
  signals = generate_signals(T)
  noise = generate_noise(T)
  
  X = cbind(signals, noise)
  colnames(X) <- paste0("V", 1:10)
  
  return(X)
}

# Function for standardizing the data.
standardize_data = function(X){
  # Center data.
  Xc = scale(X, center=TRUE, scale=FALSE)
  
  # Compute lag-zero covariance matrix.
  S0 = cov(X)

  # Perform eigen decomposition.
  eigen = eigen(S0)
  Q = eigen$vectors
  L = diag(eigen$values)
  L_inv_sqrt = sqrt(solve(L))
  
  # Compute ZCA whitening matrix.
  W = Q %*% L_inv_sqrt %*% t(Q)
  
  # Numerical verification
  check = W %*% S0 %*% t(W)
  if ((all.equal(check, diag(10), tolerance=1e-10)) == FALSE){
    print("Warning! Significantly high magnitude off-diagonals found while whitening.")
  }
  
  # Compute standardized data.
  X_st = Xc %*% W
  
  return(list(X_st=X_st, W=W))
}

# Function for computing lagged covariance matrices. Force symmetry for only real eigenvalues.
lag_cov = function(X_st, tau){
  N = nrow(X_st)
  C = cov(X_st[1:(N-tau),], X_st[(1+tau):N,])
  C_sym = 1/2 * (C + t(C))
  return(C_sym)
}

# AMUSE implementation.
AMUSE = function(SIGMA_tau){
  # Determine the amount of different lags.
  n_tau = length(SIGMA_tau)
  tau = 1:n_tau
  
  # Eigen decomposition for lagged covariance matrices.
  eigen_tau = lapply(SIGMA_tau, function(SIGMA_tau) {eigen(SIGMA_tau, symmetric=TRUE)})

  # Diagonalize lagged covariance matrices to find U  
  U = lapply(eigen_tau, function(eigen_tau) {eigen_tau$vectors})
  LAMBDA = lapply(eigen_tau, function(eigen_tau) {eigen_tau$values})[[1]]
  LAMBDA = rev(sort(LAMBDA**2))

  return(list(U=U[[1]], LAMBDA=LAMBDA))
}

SOBI = function(SIGMA_tau) {
  # Produce a S_tau 3D array for ease of computation in JADE
  SIGMA_tau_array = simplify2array(SIGMA_tau)
 
  # Use JADEs rjd function to jointly diagonalize the lagged covariances
  result = JADE::frjd(SIGMA_tau_array, maxiter = 1000)
  U = result$V
  
  D = lapply(tau, function(i) {t(U) %*% S_tau[[i]] %*% U})
  
  LAMBDA = Reduce('+', lapply(D, function(D) {diag(D)**2}))
  LAMBDA = rev(sort(LAMBDA))

  return(list(U=U, LAMBDA=LAMBDA))
}

# Function for generating fixed block bootstrap sample. X is the data and L is
# the fixed block length.
fixed_block_bootstrap = function(X, L) {
  # Number of observations.
  N = nrow(X)
  
  # Number of blocks.
  n_blocks = ceiling(N / L)
  
  # Random starting points.
  starts = sample.int(N - L + 1, n_blocks, replace = TRUE)
  
  # Construct randomly ordered indices.
  idx = unlist(lapply(starts, function(starts) {starts:(starts + L - 1)}))

  return(X[idx[seq_len(N)], , drop=FALSE])
}

# Function for generating a stationary bootstrap sample. X is the data and LE is
# the expected block length.
stationary_bootstrap = function(X, LE) {
  # Number of observations.
  N = nrow(X) 
  
  # Empty array for the soon to be generated indexes.
  idx = integer(N)
  idx[1] = sample.int(N, 1) # Random starting position.
  
  for (i in 2:N){
    if (runif(1) < 1 / LE){
      idx[i] = sample.int(N, 1) # Start new block.
    }
    else {
      idx[i] = (idx[i-1] %% N) + 1 # Continue. Go back to X_1 if N is exceeded.
    }
  }
  
  return(X[idx[seq_len(N)], , drop=FALSE])
}

# Function for computing the f0 ladle quantity.
compute_f02 = function(U, U_boot){
  f0 = sapply(1:10, function(k) {
    Bk = U[, 1:k, drop = FALSE]
    
    mean(sapply(U_boot, function(Ub) {
      Bk_boot = Ub[, 1:k, drop = FALSE]
      1 - abs(det(crossprod(Bk, Bk_boot)))
    }))
  })
  
  return(f0)
}

compute_f0 <- function(U, U_boot) {
  p <- ncol(U)
  
  f0 <- c(0, sapply(1:(p-1), function(k) {
    Bk <- U[, 1:k, drop = FALSE]
    
    mean(sapply(U_boot, function(Ub) {
      Bk_boot <- Ub[, 1:k, drop = FALSE]
      
      1 - abs(det(t(Bk) %*% Bk_boot))
    }))
  }))
  
  return(f0)
}

# Orthonormalize unmixing matrices and calculate crossproduct determinant to get
# an estimation of ladle quantity f.
subspace_distance = function(G1, G2, k) {
  
  Q1 = qr.Q(qr(G1[, 1:k, drop = FALSE]))
  Q2 = qr.Q(qr(G2[, 1:k, drop = FALSE]))
  
  return(1 - abs(det(t(Q1) %*% Q2)))
}

repetitions = 10
tau_AMUSE = c(1)
T = 1000
B = 200
L = 20

phi_AMUSE = replicate(repetitions, numeric(10), simplify = FALSE)
f_AMUSE_fixed = replicate(repetitions, numeric(10), simplify = FALSE)
g_AMUSE_fixed = replicate(repetitions, numeric(10), simplify = FALSE)

phi_SOBI = replicate(repetitions, numeric(10), simplify = FALSE)
f_SOBI_fixed = replicate(repetitions, numeric(10), simplify = FALSE)
g_SOBI_fixed = replicate(repetitions, numeric(10), simplify = FALSE)

f_AMUSE_stationary = replicate(repetitions, numeric(10), simplify = FALSE)
g_AMUSE_stationary = replicate(repetitions, numeric(10), simplify = FALSE)

f_SOBI_stationary = replicate(repetitions, numeric(10), simplify = FALSE)
g_SOBI_stationary = replicate(repetitions, numeric(10), simplify = FALSE)

for (rep in 1:repetitions){

  # Generate data.
  X = generate_data(T)

  # Standardize data.
  result = standardize_data(X)
  X_st = result$X_st
  W = result$W

  # Compute lagged covariance matrices at 20 different lags.
  tau = 1:12
  S_tau = lapply(tau, function(tau){
    lag_cov(X_st, tau)
  })
  
  # Apply AMUSE.
  AMUSE_result = AMUSE(S_tau[tau_AMUSE])
  U_AMUSE = AMUSE_result$U
  LAMBDA_AMUSE = AMUSE_result$LAMBDA
  
  # Apply SOBI.
  SOBI_result = SOBI(S_tau)
  U_SOBI = SOBI_result$U
  LAMBDA_SOBI = SOBI_result$LAMBDA
  
  # Compute ladle quantity phi for both methods.
  phi_AMUSE[[rep]] = LAMBDA_AMUSE / (1 + sum(LAMBDA_AMUSE))
  phi_SOBI[[rep]] = LAMBDA_SOBI / (1 + sum(LAMBDA_SOBI))
  
  # Initialize empty arrays to be populated with bootstrap sample results.
  U_AMUSE_boot = replicate(B, matrix(0, 10, 10), simplify = FALSE)
  U_SOBI_boot = replicate(B, matrix(0, 10, 10), simplify = FALSE)
  W_boot = replicate(B, matrix(0, 10, 10), simplify = FALSE)

  # Fixed block bootstrap.
  for (b in 1:B){
    X_boot = fixed_block_bootstrap(X, L)
    result = standardize_data(X_boot)
    X_boot_st = result$X_st
    W_boot[[b]] = result$W
    
    S_boot_tau = lapply(tau, function(tau){
      lag_cov(X_boot_st, tau)
    })
    
    AMUSE_boot = AMUSE(S_boot_tau[tau_AMUSE])
    SOBI_boot = SOBI(S_boot_tau)
    
    U_AMUSE_boot[[b]] = AMUSE_boot$U
    U_SOBI_boot[[b]] = SOBI_boot$U
  }
  
  Gamma_AMUSE = W %*% U_AMUSE
  Gamma_SOBI = W %*% U_SOBI
  
  Gamma_AMUSE_boot = Map(function(W, U) {W %*% U}, W_boot, U_AMUSE_boot)
  Gamma_SOBI_boot = Map(function(W, U) {W %*% U}, W_boot, U_SOBI_boot)
  
  # Compute ladle quantity f for both AMUSE and SOBI.
  #f0_AMUSE = compute_f0(U_AMUSE, U_AMUSE_boot)
  p = 10
  k = 0:9
  
  f0_AMUSE = sapply(k, function(k) {
    mean(sapply(Gamma_AMUSE_boot, function(Gb) {
    subspace_distance(Gamma_AMUSE, Gb, k)
    }))})
  
  f_AMUSE_fixed[[rep]] = f0_AMUSE / (1 + sum(f0_AMUSE))
  
  f0_SOBI = compute_f0(U_SOBI, U_SOBI_boot)
  
  #f0_SOBI = sapply(k, function(k) {
  #  mean(sapply(Gamma_SOBI_boot, function(Gb) {
  #  subspace_distance(Gamma_SOBI, Gb, k)
  #  }))
  #})
  
  f_SOBI_fixed[[rep]] = f0_SOBI / (1 + sum(f0_SOBI))
  
  # Compute ladle quantity g for both AMUSE and SOBI.
  g_AMUSE_fixed[[rep]] = f_AMUSE_fixed[[rep]] + phi_AMUSE[[rep]]
  g_SOBI_fixed[[rep]] = f_SOBI_fixed[[rep]] + phi_SOBI[[rep]]
  
  # Initialize empty arrays to be populated with bootstrap sample results.
  U_AMUSE_boot = replicate(B, matrix(0, 10, 10), simplify = FALSE)
  U_SOBI_boot = replicate(B, matrix(0, 10, 10), simplify = FALSE)
  W_boot = replicate(B, matrix(0, 10, 10), simplify = FALSE)
  
  # Stationary bootstrap.
  for (b in 1:B){
    X_boot = stationary_bootstrap(X , L)
    result = standardize_data(X_boot)
    X_boot_st = result$X_st
    W_boot[[b]] = result$W
    
    S_boot_tau = lapply(tau, function(tau){
      lag_cov(X_boot_st, tau)
    })
    
    AMUSE_boot = AMUSE(S_boot_tau[tau_AMUSE])
    SOBI_boot = SOBI(S_boot_tau)
    
    U_AMUSE_boot[[b]] = AMUSE_boot$U
    U_SOBI_boot[[b]] = SOBI_boot$U
  }
  
  Gamma_AMUSE_boot = Map(function(W, U) {W %*% U}, W_boot, U_AMUSE_boot)
  
  Gamma_SOBI_boot = Map(function(W, U) {W %*% U}, W_boot, U_SOBI_boot)
  
  # Compute ladle quantity f for both AMUSE and SOBI.
  p = 10
  k = 0:9
  
  #f0_AMUSE = compute_f0(U_AMUSE, U_AMUSE_boot)
  f0_AMUSE = sapply(k, function(k) {
    mean(sapply(Gamma_AMUSE_boot, function(Gb) {
    subspace_distance(Gamma_AMUSE, Gb, k)
    }))
  })
  f_AMUSE_stationary[[rep]] = f0_AMUSE / (1 + sum(f0_AMUSE))
  
  f0_SOBI = compute_f0(U_SOBI, U_SOBI_boot)
  #f0_SOBI = sapply(k, function(k) {
  #  mean(sapply(Gamma_SOBI_boot, function(Gb) {
  #  subspace_distance(t(Gamma_SOBI), Gb, k)
  #  }))
  #})
  f_SOBI_stationary[[rep]] = f0_SOBI / (1 + sum(f0_SOBI))

  # Compute ladle quantity g for both AMUSE and SOBI.
  g_AMUSE_stationary[[rep]] = f_AMUSE_stationary[[rep]] + phi_AMUSE[[rep]]
  g_SOBI_stationary[[rep]] = f_SOBI_stationary[[rep]] + phi_SOBI[[rep]]
  
}

plot(0:9, Reduce("+", g_SOBI_stationary)/repetitions, type = "l", ylim=c(0,0.7), col="green")
lines(0:9, Reduce("+", phi_SOBI)/repetitions, type="l", col="red")
lines(0:9, Reduce("+", f_SOBI_stationary)/repetitions, type="l", col="black")




