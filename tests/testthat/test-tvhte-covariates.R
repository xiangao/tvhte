# Smoke test for Phase 2: covariates X_it'beta added to the outcome.
# Two covariates with known coefficients; check the QMLE recovers beta
# and that the dynamic parameters still come back close to truth.

test_that("tvhte with covariates recovers beta and dynamic params", {
  true_beta <- c(0.7, -0.4)
  sim <- simulate_tvhte(N = 1000, T = 6, t0 = 3, J = 3,
                        rho_Y = 0.5, rho_delta = 0.7,
                        sigma_U = 1, sigma_eps = 0.3,
                        mu_alpha = 0, mu_delta0 = 1,
                        sigma_alpha = 0.5, sigma_delta0 = 0.4,
                        cor_alpha_delta = 0.2,
                        beta = true_beta, seed = 3)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X)

  expect_s3_class(fit, "tvhte")
  expect_equal(fit$convergence, 0)
  expect_equal(length(fit$beta), 2)
  expect_lt(max(abs(fit$beta - true_beta)), 0.10)

  # Dynamics still ok
  expect_lt(abs(fit$theta$rho_Y     - 0.5), 0.10)
  expect_lt(abs(fit$theta$rho_delta - 0.7), 0.15)
})

test_that("tvhte without X reduces to Phase 1 behaviour", {
  sim <- simulate_tvhte(N = 600, T = 6, t0 = 3, J = 2, seed = 4)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J)
  expect_equal(length(fit$beta), 0)
  expect_equal(fit$K, 0)
})
