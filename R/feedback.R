#' Estimate the homogeneous covariate feedback process
#'
#' Implements the covariate-feedback estimation step of Botosaru and
#' Liu (2026): under the homogeneous-feedback assumption that
#' \deqn{f_t(X_{it} | I_i^{t-1}, \lambda_i) = f_t(X_{it} | I_i^{t-1})}
#' (the covariate adjustment rule is the same across units conditional
#' on observable history), the feedback process can be estimated
#' separately from the structural outcome model. The likelihood factors
#' into a structural piece (for `Y | X, history`) and a feedback piece
#' (for `X | history`); this function fits the latter.
#'
#' Phase 5 scope: a single covariate (`K = 1`), linear AR(1) feedback
#' with one lag of `Y` and one lag of `X`:
#' \deqn{X_{it} = \gamma_0 + \gamma_Y Y_{i,t-1} + \gamma_X X_{i,t-1} + \eta_{it}.}
#' Fit by pooled OLS across unit-time observations with `t >= 2`.
#'
#' @param Y N x T outcome matrix (post-baseline).
#' @param Y0 Length-N baseline outcome.
#' @param X N x T x 1 covariate array. (Multi-covariate feedback is a
#'   future extension; pass a single covariate for now.)
#' @param X0 Length-N baseline covariate value, used as the lag at t = 1.
#'
#' @return A list of class `"tvhte_feedback"` with:
#' \describe{
#'   \item{coef}{Named numeric: `intercept, gamma_Y, gamma_X`.}
#'   \item{sigma_eta}{Residual SD of the feedback equation.}
#'   \item{fit}{The underlying `lm` object.}
#' }
#'
#' @references
#' Botosaru, Irene and Laura Liu (2026). "Event Studies with Feedback."
#' AEA Papers and Proceedings 116: 70-74.
#'
#' @examples
#' # Simulate from a model with covariate feedback, then estimate it back.
#' sim <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 2,
#'                       beta = 0.4,
#'                       feedback_gamma = c(0.2, 0.3, 0.5, 0.4),
#'                       seed = 51)
#' fb <- fit_feedback(sim$Y, sim$Y0, sim$X, sim$X0)
#' print(fb)
#' fb$coef
#'
#' @export
fit_feedback <- function(Y, Y0, X, X0) {
  if (length(dim(X)) != 3 || dim(X)[3] != 1L)
    stop("fit_feedback: Phase 5 scope is K = 1 covariate")
  N <- nrow(Y); T <- ncol(Y)
  if (length(Y0) != N || length(X0) != N)
    stop("fit_feedback: Y0 / X0 lengths must equal nrow(Y)")

  # Stack (X_it, Y_{i,t-1}, X_{i,t-1}) across all i, t in 1..T
  X_flat <- as.vector(X[, , 1])                   # length N*T, row by t
  Y_lag  <- as.vector(cbind(Y0, Y[, -T]))         # N x T matrix of lags
  X_lag  <- as.vector(cbind(X0, X[, -T, 1]))

  d <- data.frame(X = X_flat, Y_lag = Y_lag, X_lag = X_lag)
  fit <- lm(X ~ Y_lag + X_lag, data = d)
  co <- coef(fit)
  sigma_eta <- sqrt(mean(residuals(fit) ^ 2))     # MLE-style RSS / n

  structure(
    list(coef = c(intercept = unname(co[1]),
                  gamma_Y   = unname(co[2]),
                  gamma_X   = unname(co[3])),
         sigma_eta = sigma_eta,
         fit = fit),
    class = "tvhte_feedback"
  )
}

#' @export
print.tvhte_feedback <- function(x, digits = 4, ...) {
  cat("Homogeneous covariate-feedback process (Botosaru-Liu 2026)\n")
  cat(sprintf("  X_it = %s + %s * Y_{i,t-1} + %s * X_{i,t-1} + eta_it\n",
              format(x$coef[1], digits = digits),
              format(x$coef[2], digits = digits),
              format(x$coef[3], digits = digits)))
  cat(sprintf("  sigma_eta = %s\n", format(x$sigma_eta, digits = digits)))
  invisible(x)
}

#' Simulate counterfactual paths under tvhte + feedback fits
#'
#' Implements Algorithm 1 of Botosaru and Liu (2026): given fits of the
#' outcome model (`tvhte`) and feedback model (`fit_feedback`), simulate
#' counterfactual joint trajectories `{Y_i^{T,*}, X_i^{T,*}}` under
#' alternative treatment timing or initial conditions.
#'
#' Latent `lambda_i = (alpha_i, delta_{i0})` is drawn from the estimated
#' Gaussian prior; event-time effects evolve via the AR(1) on
#' `delta_{ij}`; covariates evolve via the estimated feedback process;
#' outcomes accumulate the direct effect (`delta_{ij}`) plus the
#' indirect effect (`X_it' beta`).
#'
#' @param fit A `tvhte` object.
#' @param fit_fb A `tvhte_feedback` object.
#' @param t0_star Counterfactual treatment timing (scalar applied to all
#'   simulated units, or length-`N_star` vector). Use `Inf` for never-
#'   treated counterfactuals.
#' @param N_star Number of counterfactual units to simulate.
#' @param Y0_star,X0_star Counterfactual baseline values (length-`N_star`
#'   each). Default draws standard normals.
#' @param seed Optional RNG seed.
#'
#' @return A list with `Y` (`N_star x T` matrix of counterfactual
#'   outcomes), `X` (`N_star x T x 1` array), and `lambda` (`N_star x 2`
#'   matrix of drawn latent effects).
#'
#' @examples
#' # Fit the model, fit the feedback, then simulate counterfactual paths
#' # under shifted treatment timing.
#' sim <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 2, beta = 0.4,
#'                       feedback_gamma = c(0.2, 0.3, 0.5, 0.4), seed = 1)
#' fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X,
#'              compute_se = FALSE)
#' fb  <- fit_feedback(sim$Y, sim$Y0, sim$X, sim$X0)
#'
#' # Counterfactual: treatment delayed from t=3 to t=5
#' cf <- simulate_counterfactual(fit, fb, t0_star = 5,
#'                                N_star = 200, seed = 99)
#' colMeans(cf$Y)   # average counterfactual outcome path
#'
#' @export
simulate_counterfactual <- function(fit, fit_fb, t0_star,
                                    N_star = 500,
                                    Y0_star = NULL, X0_star = NULL,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!inherits(fit, "tvhte"))
    stop("simulate_counterfactual: fit must be a tvhte object")
  if (!inherits(fit_fb, "tvhte_feedback"))
    stop("simulate_counterfactual: fit_fb must be a tvhte_feedback object")
  if (length(fit$beta) != 1L)
    stop("simulate_counterfactual: Phase 5 scope is K = 1 covariate")

  T <- fit$T; J <- fit$J
  if (length(t0_star) == 1L) t0_star <- rep(t0_star, N_star)
  if (length(t0_star) != N_star)
    stop("simulate_counterfactual: t0_star length must be 1 or N_star")

  if (is.null(Y0_star)) Y0_star <- rnorm(N_star)
  if (is.null(X0_star)) X0_star <- rnorm(N_star)

  # Draw lambda from the estimated Gaussian prior
  mu <- c(fit$prior$mu_alpha, fit$prior$mu_delta0)
  Sigma <- matrix(c(fit$prior$sigma_alpha2,
                    fit$prior$cov_alpha_delta,
                    fit$prior$cov_alpha_delta,
                    fit$prior$sigma_delta0_2), 2, 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(2 * N_star), N_star, 2)
  lambda <- t(t(Z %*% L) + mu)
  alpha  <- lambda[, 1]
  delta0 <- lambda[, 2]

  # Build per-unit event-time effects via AR(1) on delta
  sigma_eps <- sqrt(fit$theta$sigma_eps2)
  delta <- matrix(NA_real_, N_star, J + 1)
  delta[, 1] <- delta0
  if (J >= 1) {
    eps <- matrix(rnorm(N_star * J, sd = sigma_eps), N_star, J)
    for (j in 1:J) delta[, j + 1] <- fit$theta$rho_delta * delta[, j] + eps[, j]
  }

  beta <- fit$beta[1]
  sigma_U <- sqrt(fit$theta$sigma_U2)
  rho_Y <- fit$theta$rho_Y
  g0 <- fit_fb$coef[1]; gY <- fit_fb$coef[2]; gX <- fit_fb$coef[3]
  sigma_eta <- fit_fb$sigma_eta

  Y <- matrix(NA_real_, N_star, T)
  X <- array(NA_real_, dim = c(N_star, T, 1))
  for (t in 1:T) {
    Y_lag <- if (t == 1) Y0_star else Y[, t - 1]
    X_lag <- if (t == 1) X0_star else X[, t - 1, 1]
    # Feedback step: draw X_t from its estimated dynamics
    X[, t, 1] <- g0 + gY * Y_lag + gX * X_lag + rnorm(N_star, sd = sigma_eta)
    # Treatment indicator at this t for each unit
    j_i <- t - t0_star
    in_window <- is.finite(j_i) & j_i >= 0 & j_i <= J
    trt <- numeric(N_star)
    if (any(in_window))
      trt[in_window] <- delta[cbind(which(in_window), j_i[in_window] + 1)]
    # Outcome: direct (alpha + trt) + indirect (X_t * beta)
    Y[, t] <- rho_Y * Y_lag + alpha + trt +
              X[, t, 1] * beta + rnorm(N_star, sd = sigma_U)
  }

  list(Y = Y, X = X, lambda = lambda, t0 = t0_star)
}
