# Smoke tests: verify the code runs end-to-end and gives output in the
# right ballpark. Use small N + compute_se = FALSE for speed. Tight
# statistical validation lives in inst/sims/mc_validation.R.

test_that("tvhte runs and returns roughly the right parameters", {
  sim <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 3,
                        rho_Y = 0.5, rho_delta = 0.7,
                        sigma_alpha = 0.5, sigma_delta0 = 0.4,
                        seed = 1)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J,
               compute_se = FALSE)
  expect_s3_class(fit, "tvhte")
  expect_equal(fit$convergence, 0)
  expect_lt(abs(fit$theta$rho_Y     - 0.5), 0.15)
  expect_lt(abs(fit$theta$rho_delta - 0.7), 0.20)
})

test_that("posterior delta trajectory follows the AR decay", {
  sim <- simulate_tvhte(N = 500, T = 7, t0 = 3, J = 4, seed = 2,
                        rho_Y = 0.4, rho_delta = 0.6,
                        sigma_alpha = 0.6, sigma_delta0 = 0.6)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J,
               compute_se = FALSE)
  # Posterior alpha_hat should correlate with truth (not high bar)
  expect_gt(cor(fit$lambda_hat[, "alpha"], sim$lambda[, 1]), 0.5)
  expect_gt(cor(fit$lambda_hat[, "delta0"], sim$lambda[, 2]), 0.3)
})
