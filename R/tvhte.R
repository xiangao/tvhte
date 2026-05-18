#' Two-step TV-HTE estimator (Botosaru and Liu 2025)
#'
#' Fits the linear-Gaussian dynamic-panel model with correlated random
#' coefficients
#' \deqn{Y_{it} = \rho_Y Y_{i,t-1} + \alpha_i + \sum_j D_{it}^j \delta_{ij} + U_{it}}
#' \deqn{\delta_{ij} = \rho_\delta \delta_{i,j-1} + \varepsilon_{ij},\; j \geq 1}
#' under common treatment timing. The latent `lambda_i = (alpha_i, delta_{i0})`
#' is integrated out under a Gaussian working assumption; in step 1 we
#' estimate the common parameters (and the lambda-prior parameters) by
#' QMLE on the resulting marginal multivariate-normal log-likelihood; in
#' step 2 we form the closed-form Gaussian posterior mean of
#' `(alpha_i, delta_{i0})` for each unit. The Botosaru-Liu (2025) paper
#' shows that the QMLE remains consistent under misspecification of the
#' Gaussian working assumption on lambda.
#'
#' @param Y N x T matrix of post-baseline outcomes.
#' @param Y0 Length-N vector of baseline outcomes `Y_{i,0}`.
#' @param t0 Treatment period (integer, in `1:T`).
#' @param J Maximum event time observed in-window (in `0:(T - t0)`).
#' @param X Optional `N x T x K` array of strictly exogenous covariates. If
#'   supplied, the model adds `X_{it}'beta` to the outcome equation and
#'   estimates `beta` jointly with the other parameters.
#' @param init Optional named list of starting values for the optimizer
#'   (overrides any of the defaults computed from the data).
#' @param control Passed to `optim`.
#'
#' @return A list of class `"tvhte"` with:
#' \describe{
#'   \item{theta}{Estimated common parameters
#'     `c(rho_Y, rho_delta, sigma_U, sigma_eps)`.}
#'   \item{prior}{Estimated Gaussian prior on `lambda_i`.}
#'   \item{loglik}{Maximised marginal log-likelihood.}
#'   \item{lambda_hat}{N x 2 matrix of posterior means of
#'     `(alpha_i, delta_{i0})`.}
#'   \item{delta_path}{N x (J+1) matrix of posterior-mean event-time
#'     trajectories. Column `j+1` is `E[delta_{i,j} | data]`.}
#'   \item{convergence}{`optim`'s return code; 0 means converged.}
#' }
#'
#' @references
#' Botosaru, Irene and Laura Liu (2025). "Time-Varying Heterogeneous
#' Treatment Effects in Event Studies." arXiv:2509.13698.
#'
#' @export
tvhte <- function(Y, Y0, t0, J, X = NULL, init = NULL,
                  control = list(maxit = 500)) {

  if (!is.matrix(Y)) stop("tvhte: Y must be a matrix (N x T)")
  N <- nrow(Y); T <- ncol(Y)
  if (length(Y0) != N) stop("tvhte: Y0 length must equal nrow(Y)")
  stopifnot(t0 >= 1, t0 <= T, J >= 0, J <= T - t0)

  K <- 0L
  if (!is.null(X)) {
    if (length(dim(X)) != 3 || dim(X)[1] != N || dim(X)[2] != T)
      stop("tvhte: X must be an N x T x K array")
    K <- dim(X)[3]
  }

  # --- parameterisation ----------------------------------------------------
  # Internal vector (real-valued): unconstrained reparam of
  #   rho_Y, rho_delta   (in (-1, 1))         via tanh
  #   sigma_U2, sigma_eps2 (> 0)              via exp
  #   sigma_alpha2, sigma_delta0_2 (> 0)      via exp
  #   cor_alpha_delta (in (-1, 1))            via tanh
  #   mu_alpha, mu_delta0  (real)             identity

  par_init <- if (!is.null(init)) init else {
    list(rho_Y = 0.3, rho_delta = 0.3,
         log_sigma_U2 = log(var(as.vector(Y))),
         log_sigma_eps2 = log(0.25 * var(as.vector(Y))),
         mu_alpha = mean(Y), mu_delta0 = 0,
         log_sigma_alpha2 = log(0.25 * var(as.vector(Y))),
         log_sigma_delta0_2 = log(0.25 * var(as.vector(Y))),
         z_cor_alpha_delta = 0)
  }
  beta_init <- if (!is.null(init$beta)) init$beta else rep(0, K)
  par_vec <- c(atanh(par_init$rho_Y),
               atanh(par_init$rho_delta),
               par_init$log_sigma_U2,
               par_init$log_sigma_eps2,
               par_init$mu_alpha,
               par_init$mu_delta0,
               par_init$log_sigma_alpha2,
               par_init$log_sigma_delta0_2,
               par_init$z_cor_alpha_delta,
               beta_init)

  unpack <- function(p) {
    list(theta = list(rho_Y = tanh(p[1]),
                      rho_delta = tanh(p[2]),
                      sigma_U2 = exp(p[3]),
                      sigma_eps2 = exp(p[4])),
         prior = list(mu_alpha = p[5],
                      mu_delta0 = p[6],
                      sigma_alpha2 = exp(p[7]),
                      sigma_delta0_2 = exp(p[8]),
                      cov_alpha_delta = tanh(p[9]) *
                        sqrt(exp(p[7]) * exp(p[8]))),
         beta = if (K > 0) p[10:(9 + K)] else numeric(0))
  }

  # --- design (same for every unit under common timing) --------------------
  build_ds <- function(theta) .unit_design(T, t0, J, theta$rho_Y, theta$rho_delta)

  # --- per-unit X %*% beta loadings, A %*% (X_i beta), recomputed only when
  #     rho_Y or beta changes within optim. We keep it inside nll for clarity.
  .unit_Xbeta <- function(ds, beta) {
    if (K == 0) return(replicate(N, NULL, simplify = FALSE))
    lapply(1:N, function(i) {
      Xi <- matrix(X[i, , ], nrow = T)         # T x K
      drop(ds$A %*% (Xi %*% beta))
    })
  }

  # --- negative log-likelihood ---------------------------------------------
  nll <- function(p) {
    par <- unpack(p)
    ds <- build_ds(par$theta)
    XB <- .unit_Xbeta(ds, par$beta)
    ll <- 0
    for (i in 1:N) {
      mom <- .unit_moments(ds, Y0[i], par$theta, par$prior, Xi_beta = XB[[i]])
      ll_i <- .dmvnorm_log(Y[i, ], mom$mean, mom$cov)
      if (is.na(ll_i) || !is.finite(ll_i)) return(.Machine$double.xmax / 2)
      ll <- ll + ll_i
    }
    -ll
  }

  fit <- optim(par_vec, nll, method = "BFGS", control = control)
  par <- unpack(fit$par)

  # --- step 2: posterior means of lambda_i ---------------------------------
  ds <- build_ds(par$theta)
  XB <- .unit_Xbeta(ds, par$beta)
  Sigma_lambda <- matrix(c(par$prior$sigma_alpha2,
                           par$prior$cov_alpha_delta,
                           par$prior$cov_alpha_delta,
                           par$prior$sigma_delta0_2), 2, 2)
  mu_lambda <- c(par$prior$mu_alpha, par$prior$mu_delta0)

  # Cross-covariance between (alpha_i, delta_{i0}) and Y_i (constant across i)
  Sigma_lY <- rbind(
    ds$A1  * par$prior$sigma_alpha2 + ds$AMc * par$prior$cov_alpha_delta,
    ds$A1  * par$prior$cov_alpha_delta + ds$AMc * par$prior$sigma_delta0_2
  )

  lambda_hat <- matrix(NA_real_, N, 2)
  for (i in 1:N) {
    mom <- .unit_moments(ds, Y0[i], par$theta, par$prior, Xi_beta = XB[[i]])
    inv_S_Y <- solve(mom$cov)
    lambda_hat[i, ] <- mu_lambda + Sigma_lY %*% inv_S_Y %*% (Y[i, ] - mom$mean)
  }
  colnames(lambda_hat) <- c("alpha", "delta0")

  # Posterior event-time trajectory: E[delta_{i,j} | data] for j = 0..J.
  # Under Gaussian linear model with mean function depending linearly on
  # (alpha, delta_0, eps_1, ..., eps_J), the conditional means of the
  # innovations eps_{i,j} are also linear in (Y - E[Y]). Build cross-cov
  # of the full latent vector and project.
  delta_path <- matrix(NA_real_, N, J + 1)
  delta_path[, 1] <- lambda_hat[, 2]   # E[delta_{i0} | data]
  if (J >= 1) {
    # Cross-cov of eps_i (J x 1) with Y_i (T x 1): eps enters Y via AMW with
    # variance sigma_eps2 * I_J, independent of alpha and delta0.
    Sigma_epsY <- t(ds$AMW) * par$theta$sigma_eps2     # J x T
    eps_hat <- matrix(NA_real_, N, J)
    for (i in 1:N) {
      mom <- .unit_moments(ds, Y0[i], par$theta, par$prior, Xi_beta = XB[[i]])
      inv_S_Y <- solve(mom$cov)
      eps_hat[i, ] <- Sigma_epsY %*% inv_S_Y %*% (Y[i, ] - mom$mean)
    }
    # delta_{i,j} = rho_delta^j * delta_{i0} + sum_{k=1}^j rho_delta^{j-k} eps_{i,k}
    rho_d <- par$theta$rho_delta
    for (j in 1:J) {
      delta_path[, j + 1] <- rho_d^j * lambda_hat[, 2] +
        as.numeric(eps_hat[, 1:j, drop = FALSE] %*% rho_d^((j - 1):0))
    }
  }
  colnames(delta_path) <- paste0("delta_", 0:J)

  structure(
    list(theta = par$theta, prior = par$prior, beta = par$beta,
         loglik = -fit$value, convergence = fit$convergence,
         lambda_hat = lambda_hat, delta_path = delta_path,
         t0 = t0, J = J, N = N, T = T, K = K,
         call = match.call()),
    class = "tvhte"
  )
}

#' @export
print.tvhte <- function(x, digits = 4, ...) {
  cat("Time-Varying Heterogeneous Treatment Effects (Botosaru-Liu 2025)\n")
  cat(sprintf("  N = %d units; T = %d periods; t0 = %d; J = %d\n",
              x$N, x$T, x$t0, x$J))
  cat(sprintf("  log-likelihood = %.3f   convergence = %d\n\n",
              x$loglik, x$convergence))
  cat("Common parameters (theta):\n")
  cat(sprintf("  rho_Y      = %s\n", format(x$theta$rho_Y, digits = digits)))
  cat(sprintf("  rho_delta  = %s\n", format(x$theta$rho_delta, digits = digits)))
  cat(sprintf("  sigma_U    = %s\n", format(sqrt(x$theta$sigma_U2), digits = digits)))
  cat(sprintf("  sigma_eps  = %s\n", format(sqrt(x$theta$sigma_eps2), digits = digits)))
  cat("Prior on lambda_i = (alpha, delta_0):\n")
  cat(sprintf("  mu = (%s, %s)\n",
              format(x$prior$mu_alpha, digits = digits),
              format(x$prior$mu_delta0, digits = digits)))
  cat(sprintf("  sd = (%s, %s),  cor = %s\n",
              format(sqrt(x$prior$sigma_alpha2), digits = digits),
              format(sqrt(x$prior$sigma_delta0_2), digits = digits),
              format(x$prior$cov_alpha_delta /
                       sqrt(x$prior$sigma_alpha2 * x$prior$sigma_delta0_2),
                     digits = digits)))
  if (length(x$beta) > 0) {
    cat(sprintf("\nCovariate coefficients (beta):\n"))
    cat(sprintf("  %s\n",
                paste0("beta[", seq_along(x$beta), "] = ",
                       format(x$beta, digits = digits),
                       collapse = "   ")))
  }
  cat(sprintf("\nMean posterior event-time effects (across units):\n"))
  cat(sprintf("  %s\n", paste(colnames(x$delta_path),
                              format(colMeans(x$delta_path), digits = digits),
                              sep = " = ", collapse = "   ")))
  invisible(x)
}
