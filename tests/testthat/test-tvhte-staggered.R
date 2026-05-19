# Phase 3: staggered adoption.

test_that("staggered cohorts including Inf are handled", {
  N <- 400
  set.seed(20)
  cohorts <- sample(c(3, 4, 5, Inf), size = N, replace = TRUE,
                    prob = c(0.3, 0.25, 0.2, 0.25))
  sim <- simulate_tvhte(N = N, T = 7, t0 = cohorts, J = 2, seed = 21)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J,
               compute_se = FALSE)
  expect_equal(fit$convergence, 0)
  expect_equal(length(fit$cohort_counts), 4)
  expect_true("Inf" %in% names(fit$cohort_counts))
  expect_lt(abs(fit$theta$rho_Y     - 0.5), 0.20)
})

test_that("scalar t0 still works (backwards compatibility)", {
  sim <- simulate_tvhte(N = 300, T = 5, t0 = 3, J = 2, seed = 22)
  fit <- tvhte(sim$Y, sim$Y0, t0 = 3, J = sim$J,
               compute_se = FALSE)
  expect_equal(length(unique(fit$t0)), 1)
})

test_that("staggered with covariates works", {
  N <- 500
  set.seed(30)
  cohorts <- sample(c(3, 5, Inf), size = N, replace = TRUE,
                    prob = c(0.4, 0.3, 0.3))
  true_beta <- c(0.5, -0.3)
  sim <- simulate_tvhte(N = N, T = 6, t0 = cohorts, J = 2,
                        beta = true_beta, seed = 31)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X,
               compute_se = FALSE)
  expect_lt(max(abs(fit$beta - true_beta)), 0.15)
})
