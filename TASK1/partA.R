# Read in the data (using the actual data now in the beginning just to get a better idea of what is actually supposed to happen)
# Will switch to generated data at some point
X = read.table("/home/elmer/R/computational_statistics/group_work/TASK1/time_series_data.txt")

# compute the amount of observations and variables
N = nrow(X)
p = ncol(X)

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

#
tau = 1:20
S_tau = lapply(tau, function(tau){
  lag_cov(X_st, tau)
})

eigen(S_tau[[5]])