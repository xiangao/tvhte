# Phase 3: staggered adoption.
#
# Three treated cohorts (t0 = 3, 4, 5) plus a never-treated group (t0 = Inf).
# Same true (rho_Y, rho_delta, beta) as Phase 2; check QMLE still recovers.

test_that("tvhte handles staggered cohorts + never-treated", {
  N <- 1500
  set.seed(20)
  cohorts <- sample(c(3, 4, 5, Inf), size = N, replace = TRUE,
                    prob = c(0.3, 0.25, 0.2, 0.25))
  sim <- simulate_tvhte(N = N, T = 7, t0 = cohorts, J = 2,
                        rho_Y = 0.4, rho_delta = 0.6,
                        sigma_U = 1, sigma_eps = 0.3,
                        mu_alpha = 0, mu_delta0 = 1,
                        sigma_alpha = 0.5, sigma_delta0 = 0.4,
                        cor_alpha_delta = 0.2,
                        seed = 21)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J)

  expect_s3_class(fit, "tvhte")
  expect_equal(fit$convergence, 0)
  expect_equal(length(fit$cohort_counts), 4)
  # Never-treated cohort should appear as "Inf"
  expect_true("Inf" %in% names(fit$cohort_counts))

  expect_lt(abs(fit$theta$rho_Y     - 0.4), 0.10)
  expect_lt(abs(fit$theta$rho_delta - 0.6), 0.15)
})

test_that("scalar t0 still works (backwards compatibility)", {
  sim <- simulate_tvhte(N = 400, T = 5, t0 = 3, J = 2, seed = 22)
  fit <- tvhte(sim$Y, sim$Y0, t0 = 3, J = sim$J)
  expect_s3_class(fit, "tvhte")
  expect_equal(length(unique(fit$t0)), 1)
  expect_equal(fit$t0[1], 3)
})

test_that("staggered with covariates works", {
  N <- 1000
  set.seed(30)
  cohorts <- sample(c(3, 5, Inf), size = N, replace = TRUE,
                    prob = c(0.4, 0.3, 0.3))
  true_beta <- c(0.5, -0.3)
  sim <- simulate_tvhte(N = N, T = 6, t0 = cohorts, J = 2,
                        beta = true_beta, seed = 31)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X)
  expect_lt(max(abs(fit$beta - true_beta)), 0.10)
})
