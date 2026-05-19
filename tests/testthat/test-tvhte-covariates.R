# Covariates smoke tests: verify X_it'beta path runs and beta lands
# in the right neighbourhood.

test_that("tvhte with covariates returns beta in the right ballpark", {
  true_beta <- c(0.7, -0.4)
  sim <- simulate_tvhte(N = 400, T = 6, t0 = 3, J = 3,
                        beta = true_beta, seed = 3)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X,
               compute_se = FALSE)
  expect_equal(length(fit$beta), 2)
  expect_lt(max(abs(fit$beta - true_beta)), 0.15)
})

test_that("tvhte without X reports K = 0 and empty beta", {
  sim <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 2, seed = 4)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J,
               compute_se = FALSE)
  expect_equal(length(fit$beta), 0)
  expect_equal(fit$K, 0)
})
