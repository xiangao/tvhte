# ============================================================================
# Monte Carlo validation for tvhte
# ============================================================================
#
# Verifies bias, SE accuracy, and 95% CI coverage of the QMLE step
# (rho_Y, rho_delta, beta) on the linear-Gaussian DGP across a range of
# sample sizes and treatment-effect persistence values. Also tests
# robustness when the latent lambda is non-Gaussian (the Botosaru-Liu
# QMLE is supposed to stay consistent under such misspecification).
#
# Run from the package root:
#   Rscript inst/sims/mc_validation.R [reps]
# Default reps = 100.
# ============================================================================

suppressMessages({ library(stats); library(parallel) })
source("R/design.R")
source("R/simulate.R")
source("R/tvhte.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Use most cores but leave one free for the OS
NCORES <- max(1L, detectCores() - 1L)
cat(sprintf("Using %d cores in parallel\n", NCORES))

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 100L

# ---- DGP 1: vanilla linear-Gaussian model -------------------------------
sim_linear <- function(N = 1000, seed = 1, rho_delta = 0.6) {
  simulate_tvhte(N = N, T = 6, t0 = 3, J = 3,
                 rho_Y = 0.5, rho_delta = rho_delta,
                 sigma_U = 1, sigma_eps = 0.3,
                 mu_alpha = 0, mu_delta0 = 1,
                 sigma_alpha = 0.5, sigma_delta0 = 0.4,
                 cor_alpha_delta = 0.2,
                 beta = c(0.5, -0.3), seed = seed)
}

# ---- DGP 2: bimodal lambda (mixture of two Gaussians) -------------------
# The QMLE should remain consistent for the common parameters under this
# misspecification; prior parameters are not interpretable but the
# dynamics (rho_Y, rho_delta) and beta should still recover correctly.
sim_bimodal <- function(N = 1000, seed = 1) {
  set.seed(seed)
  N <- as.integer(N)
  # Two latent groups with different means; mixture is non-Gaussian.
  z <- sample(c(-1, 1), N, replace = TRUE)
  mu_alpha_i  <- 0.7 * z
  mu_delta0_i <- 1.0 + 0.5 * z
  alpha  <- mu_alpha_i  + rnorm(N, sd = 0.3)
  delta0 <- mu_delta0_i + rnorm(N, sd = 0.3)
  T <- 6; t0 <- 3; J <- 3
  rho_Y <- 0.5; rho_delta <- 0.6
  sigma_U <- 1; sigma_eps <- 0.3
  beta <- c(0.5, -0.3)

  eps <- matrix(rnorm(N * J, sd = sigma_eps), N, J)
  delta <- matrix(NA_real_, N, J + 1); delta[, 1] <- delta0
  for (j in 1:J) delta[, j + 1] <- rho_delta * delta[, j] + eps[, j]

  Y0 <- rnorm(N)
  X  <- array(rnorm(N * T * 2), dim = c(N, T, 2))
  Y  <- matrix(NA_real_, N, T)
  for (t in 1:T) {
    Y_lag <- if (t == 1) Y0 else Y[, t - 1]
    trt <- if (t >= t0 && t - t0 <= J) delta[, t - t0 + 1] else 0
    x_eff <- as.numeric(matrix(X[, t, ], nrow = N) %*% beta)
    Y[, t] <- rho_Y * Y_lag + alpha + trt + x_eff + rnorm(N, sd = sigma_U)
  }
  list(Y = Y, Y0 = Y0, X = X, t0 = t0, J = J,
       params = list(rho_Y = rho_Y, rho_delta = rho_delta,
                     sigma_U = sigma_U, sigma_eps = sigma_eps,
                     beta = beta))
}

# ---- MC harness ---------------------------------------------------------
# Each rep is one fit; parallelise over reps using mclapply (forking, so
# the worker inherits our sourced functions automatically).
fit_one <- function(r, sim_fn, truth) {
  sim <- sim_fn(seed = r)
  fit <- tryCatch(
    tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X),
    error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0)
    return(setNames(rep(list(c(NA, NA, NA)), length(truth)), names(truth)))
  out <- vector("list", length(truth))
  names(out) <- names(truth)
  for (nm in names(truth)) {
    val <- switch(nm,
      rho_Y     = fit$theta$rho_Y,
      rho_delta = fit$theta$rho_delta,
      beta1     = fit$beta[1],
      beta2     = fit$beta[2])
    se <- if (is.null(fit$se)) NA_real_ else switch(nm,
      rho_Y     = fit$se$rho_Y     %||% NA_real_,
      rho_delta = fit$se$rho_delta %||% NA_real_,
      beta1     = if (length(fit$se$beta) >= 1) fit$se$beta[1] else NA_real_,
      beta2     = if (length(fit$se$beta) >= 2) fit$se$beta[2] else NA_real_)
    cov_hit <- if (is.na(se)) NA
               else truth[[nm]] >= val - 1.96 * se &&
                    truth[[nm]] <= val + 1.96 * se
    out[[nm]] <- c(val, se, cov_hit)
  }
  out
}

mc_run <- function(sim_fn, reps, truth) {
  res <- mclapply(seq_len(reps),
                  function(r) fit_one(r, sim_fn, truth),
                  mc.cores = NCORES, mc.preschedule = FALSE)
  out <- data.frame(param = names(truth),
                    truth = unlist(truth),
                    bias  = NA_real_, emp_sd = NA_real_,
                    mean_se = NA_real_, coverage = NA_real_)
  rownames(out) <- NULL
  for (i in seq_along(truth)) {
    nm <- names(truth)[i]
    triples <- do.call(rbind, lapply(res, `[[`, nm))
    out$bias[i]     <- mean(triples[, 1], na.rm = TRUE) - truth[[nm]]
    out$emp_sd[i]   <- sd(triples[, 1], na.rm = TRUE)
    out$mean_se[i]  <- mean(triples[, 2], na.rm = TRUE)
    out$coverage[i] <- mean(triples[, 3], na.rm = TRUE)
  }
  out
}

cat(sprintf("MC validation: %d reps per scenario\n", reps))
cat("============================================================\n\n")

truth_linear <- list(rho_Y = 0.5, rho_delta = 0.6, beta1 = 0.5, beta2 = -0.3)

cat("=== DGP 1: linear-Gaussian, N = 500 ===\n")
print(mc_run(function(seed) sim_linear(N = 500, seed = seed),
             reps, truth_linear), row.names = FALSE)

cat("\n=== DGP 1: linear-Gaussian, N = 1500 ===\n")
print(mc_run(function(seed) sim_linear(N = 1500, seed = seed),
             reps, truth_linear), row.names = FALSE)

cat("\n=== DGP 2: BIMODAL lambda (QMLE robustness), N = 1500 ===\n")
print(mc_run(function(seed) sim_bimodal(N = 1500, seed = seed),
             reps, truth_linear), row.names = FALSE)

cat("\n=== DGP 3: oscillating delta (rho_delta = -0.5), N = 1500 ===\n")
truth_osc <- list(rho_Y = 0.5, rho_delta = -0.5, beta1 = 0.5, beta2 = -0.3)
print(mc_run(function(seed) sim_linear(N = 1500, seed = seed, rho_delta = -0.5),
             reps, truth_osc), row.names = FALSE)

cat("\nInterpretation:\n")
cat("- bias -> 0 as N grows; QMLE is sqrt(N) consistent.\n")
cat("- meanSE close to empSD verifies the inverse-Hessian SE is right under\n")
cat("  the Gaussian working assumption.\n")
cat("- coverage near 0.95 for the linear-Gaussian DGP (DGP 1).\n")
cat("- DGP 2 (bimodal lambda): point estimates for rho_Y, rho_delta, beta\n")
cat("  should remain ~unbiased even though the Gaussian working assumption\n")
cat("  is wrong. Coverage may degrade since the SE is naive (no sandwich).\n")
cat("- DGP 3: oscillating delta path; estimator should recover rho_delta < 0.\n")
