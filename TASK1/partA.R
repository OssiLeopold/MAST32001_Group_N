
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
  s1 = generate_arma(ar = 0.6, ma = c(0, -0.4))
  s2 = generate_arma(ar = c(0.4, 0.1, 0.3), ma = c(0.2, 0.4))
  s3 = generate_arma(ar = 0.4)
  s4 = generate_arma(ma = c(0.5, 0.3, 0.1, 0.2))
  s5 = generate_arma(ar = c(0.1, 0.3))
  
  signals = cbind(s1, s2, s3, s4, s5)
  
  return(signals)
}

# Function to generate noise according to setting 1: Gaussian white noise.
generate_noise = function(T) {
  noise = matrix(rnorm(T * 5), nrow = T, ncol = 5)
  
  return(noise)
}

# Function for generating random observations for testing.
generate_data = function(T, setting = 1) {
  signals = generate_signals(T)
  noise = generate_noise(T)
  
  X = cbind(signals, noise)
  
  return(X)
}

# Function for standardizing the data.
standardize_data = function(X){
  # Center data.
  Xc = scale(X, center=TRUE, scale=FALSE)
  
  # Compute lag-zero covariance matrix.
  S0 = cov(Xc)
  
  # Perform eigen decomposition.
  eigen = eigen(S0)
  Q = eigen$vectors
  L = eigen$values
  L_inv_sqrt = sqrt(solve(L))
  
  # Compute ZCA whitening matrix.
  W = Q %*% L_inv_sqrt %*% t(L_inv_sqrt)
  
  # Numerical verification
  check = W %*% S0 %*% t(W)
  if (round(check, 1e-10) != diag(ncol(X))){
    print("Warning! Significantly high magnitude off-diagonals found while whitening.")
  }
  
  # Compute standardized data.
  X_st = Xc %*% W
  
  return(X_st)
}

# center the data such that each column has zero mean
Xc = scale(X, center=TRUE, scale=FALSE)

# compute the lag-zero covariance matrix
S = cov(Xc)

# perform eigen decomposition
eigen = eigen(S)
Q = eigen$vectors
QT = t(Q)
lambda = diag(eigen$values)
lambda_inv = solve(lambda)
lambda_inv_sqrt = sqrt(lambda_inv)

# compute the ZCA whitening matrix as per Nordhausen 2018
S_inv_sqrt = Q %*% lambda_inv_sqrt %*% QT

# verify that the zero-lag covariance matrix becomes the identity after whitening
S_ZCA = S_inv_sqrt %*% S %*% t(S_inv_sqrt)

# whiten the observaitons to construct x_st as per Nordhausen 2018
X_st = Xc %*% S_inv_sqrt

# define function for computing the lagged covariance matrices. Force symmetry to get only real eigenvalues.
lag_cov = function(X_st, tau){
  C = cov(X_st[1:(N-tau),], X_st[(1+tau):N,])
  C_sym = 1/2 * (C + t(C))
  return(C_sym)
}

# Compute lagged covariance matrices at 20 different lags
tau = 1:20
S_tau = lapply(tau, function(tau){
  lag_cov(X_st, tau)
})

# use AMUSE first
eigen_array = lapply(S_tau, function(S_tau){
  eigen(S_tau, symmetric=TRUE)
})

# compute the orthogonal matrix U for a range of lags
U_AMUSE_array = lapply(eigen_array, function(eigen_array){
  eigen_array$vectors
})

# compute U for range of taus
D_AMUSE = lapply(tau, function(i) {
  t(U_AMUSE_array[[i]]) %*% S_tau[[i]] %*% U_AMUSE_array[[i]]
})

L_AMUSE_array = lapply(tau, function(i){
  diag(t(U_AMUSE_array[[i]]) %*% S_tau[[i]] %*% U_AMUSE_array[[i]])**2
})

L_AMUSE_array

# Use SOBI
library(JADE)

S_tau_array = simplify2array(S_tau)

result = JADE::frjd(S_tau_array)
U_SOBI = result$V

D_SOBI = lapply(tau, function(i) {
  U_SOBI %*% S_tau[[i]] %*% t(U_SOBI)
})

L_SOBI = Reduce('+', lapply(D_SOBI, function(D){
  diag(D)**2
}))

L_SOBI

