#' Simulate from the Botosaru-Liu (2025) TV-HTE model
#'
#' Draws a dataset of N units observed at times 0, 1, ..., T according to
#' \deqn{Y_{it} = \rho_Y Y_{i,t-1} + \alpha_i + \sum_{j} D_{it}^j \delta_{ij} + U_{it}}
#' \deqn{\delta_{ij} = \rho_\delta \delta_{i,j-1} + \varepsilon_{ij}, \quad j \geq 1}
#' with `(alpha_i, delta_{i0})` jointly Gaussian.
#'
#' All units share a common treatment period `t0` and event-time window
#' `0:J` (Phase 1 scope).
#'
#' @param N Number of units.
#' @param T Number of post-baseline periods (i.e. observations at
#'   t = 1, ..., T).
#' @param t0 Treatment period (1-indexed, in `1:T`).
#' @param J Maximum event time observed in-window.
#' @param rho_Y Outcome AR(1) coefficient.
#' @param rho_delta Event-time AR(1) coefficient.
#' @param sigma_U,sigma_eps Standard deviations of `U_{it}` and `eps_{ij}`.
#' @param mu_alpha,mu_delta0,sigma_alpha,sigma_delta0,cor_alpha_delta
#'   Parameters of the Gaussian prior on `lambda_i = (alpha_i, delta_{i0})`.
#' @param Y0_mean,Y0_sd Distribution of the baseline outcome `Y_{i,0}`.
#' @param seed Optional integer seed.
#'
#' @return A list with `Y` (an N x T matrix), `Y0` (length-N baseline
#'   vector), `t0` (scalar), `J`, `lambda` (N x 2 matrix of true alpha
#'   and delta_{i0}), and `delta` (N x (J+1) matrix of true event-time
#'   effects).
#'
#' @export
simulate_tvhte <- function(N = 500, T = 6, t0 = 3, J = 3,
                           rho_Y = 0.5, rho_delta = 0.7,
                           sigma_U = 1, sigma_eps = 0.3,
                           mu_alpha = 0, mu_delta0 = 1,
                           sigma_alpha = 0.5, sigma_delta0 = 0.5,
                           cor_alpha_delta = 0,
                           Y0_mean = 0, Y0_sd = 1,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(t0 >= 1, t0 <= T, J >= 0, J <= T - t0)

  # Draw lambda_i = (alpha_i, delta_{i0}) jointly Gaussian
  Sigma_lambda <- matrix(c(sigma_alpha^2,
                           cor_alpha_delta * sigma_alpha * sigma_delta0,
                           cor_alpha_delta * sigma_alpha * sigma_delta0,
                           sigma_delta0^2), 2, 2)
  L <- chol(Sigma_lambda)
  Z <- matrix(rnorm(2 * N), N, 2)
  lambda <- t(t(Z %*% L) + c(mu_alpha, mu_delta0))
  alpha <- lambda[, 1]
  delta0 <- lambda[, 2]

  # Draw eps_{ij} for j = 1..J, then build delta_{ij}
  eps <- matrix(rnorm(N * J, sd = sigma_eps), N, J)
  delta <- matrix(NA_real_, N, J + 1)
  delta[, 1] <- delta0
  if (J >= 1) {
    for (j in 1:J) delta[, j + 1] <- rho_delta * delta[, j] + eps[, j]
  }

  # Draw Y_{i,0}
  Y0 <- rnorm(N, mean = Y0_mean, sd = Y0_sd)

  # Iterate forward
  Y <- matrix(NA_real_, N, T)
  for (t in 1:T) {
    Y_lag <- if (t == 1) Y0 else Y[, t - 1]
    trt <- if (t >= t0 && (t - t0) <= J) delta[, t - t0 + 1] else 0
    Y[, t] <- rho_Y * Y_lag + alpha + trt + rnorm(N, sd = sigma_U)
  }

  list(Y = Y, Y0 = Y0, t0 = t0, J = J,
       lambda = lambda, delta = delta,
       params = list(rho_Y = rho_Y, rho_delta = rho_delta,
                     sigma_U = sigma_U, sigma_eps = sigma_eps,
                     mu_alpha = mu_alpha, mu_delta0 = mu_delta0,
                     sigma_alpha = sigma_alpha, sigma_delta0 = sigma_delta0,
                     cor_alpha_delta = cor_alpha_delta))
}
