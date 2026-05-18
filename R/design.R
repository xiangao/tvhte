# Internal helpers building the per-unit linear-Gaussian system implied by
# Botosaru and Liu (2025), Section 2:
#
#   Y_it = rho_Y * Y_{i,t-1} + alpha_i + sum_j D_{it}^j * delta_{ij} + U_it
#   delta_{ij} = rho_delta * delta_{i,j-1} + eps_{ij},  j >= 1
#   lambda_i  = (alpha_i, delta_{i0})'
#
# Solving the AR(1) in calendar time gives
#   Y_i = rho_Y * (A r) Y_{i0}
#       + (A 1_T)        * alpha_i
#       + (A M c)        * delta_{i0}
#       + (A M W)        * eps_i
#       + A              * U_i
# where
#   A   = (I - rho_Y L)^{-1} -- T x T lower triangular, A[t,s] = rho_Y^(t-s)
#   r   = (1, 0, ..., 0)'    -- T x 1
#   M   = T x (J+1) indicator placing delta_{i,j} at calendar time t0+j
#   c   = (1, rho_delta, ..., rho_delta^J)' -- (J+1) x 1
#   W   = (J+1) x J lower triangular such that v_i = c * delta_{i0} + W eps_i
#         with W[j+1, k] = rho_delta^(j-k) for 1 <= k <= j, 0 otherwise.

# Build everything that does NOT depend on the priors, for a single unit.
.unit_design <- function(T, t0, J, rho_Y, rho_delta) {
  if (T < 1)            stop(".unit_design: T must be >= 1")
  if (J < 0)            stop(".unit_design: J must be >= 0")

  # A
  A <- matrix(0, T, T)
  pw <- rho_Y ^ (0:(T - 1))
  for (t in 1:T) for (s in 1:t) A[t, s] <- pw[t - s + 1]

  # M: place delta_{i,j} at calendar time (t0 + j); skip out-of-window cells
  M <- matrix(0, T, J + 1)
  for (j in 0:J) {
    t_idx <- t0 + j
    if (t_idx >= 1 && t_idx <= T) M[t_idx, j + 1] <- 1
  }

  c_vec <- rho_delta ^ (0:J)

  W <- matrix(0, J + 1, J)
  if (J >= 1) {
    for (j in 1:J) for (k in 1:j) W[j + 1, k] <- rho_delta ^ (j - k)
  }

  r <- c(1, rep(0, T - 1))

  AMc <- drop(A %*% M %*% c_vec)        # T-vector loading on delta_{i0}
  A1  <- drop(A %*% rep(1, T))          # T-vector loading on alpha_i
  AMW <- A %*% M %*% W                  # T x J matrix loading on eps_i

  list(A = A, M = M, c_vec = c_vec, W = W, r = r,
       AMc = AMc, A1 = A1, AMW = AMW)
}

# Per-unit marginal mean and covariance of Y_i given Y_{i0} (and X_i).
# `prior` carries the Gaussian on lambda_i. `Xi_beta`, when supplied, is
# the T-vector A %*% X_i %*% beta -- pre-multiplied so the caller can
# cache it across iterations of optim.
.unit_moments <- function(ds, Y0, theta, prior, Xi_beta = NULL) {
  T <- length(ds$A1)

  mu_Y <- theta$rho_Y * (ds$A %*% ds$r) * Y0 +
          ds$A1  * prior$mu_alpha +
          ds$AMc * prior$mu_delta0
  if (!is.null(Xi_beta)) mu_Y <- mu_Y + Xi_beta

  S_alpha <- prior$sigma_alpha2
  S_delt0 <- prior$sigma_delta0_2
  S_ad    <- prior$cov_alpha_delta
  S_eps   <- theta$sigma_eps2
  S_U     <- theta$sigma_U2

  Sigma_Y <- ds$A1 %*% t(ds$A1) * S_alpha +
             ds$AMc %*% t(ds$AMc) * S_delt0 +
             (ds$A1 %*% t(ds$AMc) + ds$AMc %*% t(ds$A1)) * S_ad +
             ds$AMW %*% t(ds$AMW) * S_eps +
             ds$A %*% t(ds$A) * S_U

  # Ridge for numerical PD-ness against degenerate parameters
  Sigma_Y <- (Sigma_Y + t(Sigma_Y)) / 2 + diag(1e-10, T)

  list(mean = as.numeric(mu_Y), cov = Sigma_Y)
}

# Multivariate normal log-density (Cholesky form, fewer ops than dmvnorm).
.dmvnorm_log <- function(y, mu, Sigma) {
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(L)) return(NA_real_)
  z <- backsolve(L, y - mu, transpose = TRUE)
  -0.5 * length(y) * log(2 * pi) - sum(log(diag(L))) - 0.5 * sum(z * z)
}
