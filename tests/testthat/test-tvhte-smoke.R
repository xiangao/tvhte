# Smoke test: simulate from the model and check the QMLE recovers the
# true common parameters within a generous tolerance, and that the
# unit-level posterior means correlate strongly with true latent effects.

test_that("tvhte QMLE recovers true common parameters", {
  sim <- simulate_tvhte(N = 800, T = 6, t0 = 3, J = 3,
                        rho_Y = 0.5, rho_delta = 0.7,
                        sigma_U = 1, sigma_eps = 0.3,
                        mu_alpha = 0, mu_delta0 = 1,
                        sigma_alpha = 0.5, sigma_delta0 = 0.4,
                        cor_alpha_delta = 0.3, seed = 1)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J)

  expect_s3_class(fit, "tvhte")
  expect_equal(fit$convergence, 0)

  # Tolerances loose enough to absorb single-replicate noise + Gaussian
  # working assumption when true lambda was Gaussian here.
  expect_lt(abs(fit$theta$rho_Y     - 0.5), 0.10)
  expect_lt(abs(fit$theta$rho_delta - 0.7), 0.15)
  expect_lt(abs(sqrt(fit$theta$sigma_U2)   - 1.0), 0.10)
})

test_that("tvhte posterior means track true latent effects", {
  sim <- simulate_tvhte(N = 1000, T = 7, t0 = 3, J = 4, seed = 2,
                        rho_Y = 0.4, rho_delta = 0.6,
                        sigma_alpha = 0.6, sigma_delta0 = 0.6)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J)

  # Posterior alpha_hat should correlate strongly with true alpha
  expect_gt(cor(fit$lambda_hat[, "alpha"], sim$lambda[, 1]), 0.7)
  # Posterior delta0_hat should correlate with true delta_{i0}
  expect_gt(cor(fit$lambda_hat[, "delta0"], sim$lambda[, 2]), 0.5)
  # Posterior delta_{i,j} trajectory should correlate with true delta column j
  for (j in 0:sim$J) {
    expect_gt(cor(fit$delta_path[, j + 1], sim$delta[, j + 1]), 0.3)
  }
})
