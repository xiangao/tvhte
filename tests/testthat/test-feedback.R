# Phase 5: homogeneous covariate feedback (Botosaru-Liu 2026).

test_that("fit_feedback recovers AR(1) feedback coefficients", {
  true_gamma <- c(intercept = 0.2, gamma_Y = 0.3,
                  gamma_X = 0.5, sigma_eta = 0.4)
  sim <- simulate_tvhte(N = 500, T = 6, t0 = 3, J = 2,
                        beta = 0.4,
                        feedback_gamma = true_gamma,
                        seed = 51)
  fb <- fit_feedback(sim$Y, sim$Y0, sim$X, sim$X0)
  expect_s3_class(fb, "tvhte_feedback")
  expect_lt(abs(fb$coef["intercept"] - true_gamma["intercept"]), 0.10)
  expect_lt(abs(fb$coef["gamma_Y"]   - true_gamma["gamma_Y"]),   0.10)
  expect_lt(abs(fb$coef["gamma_X"]   - true_gamma["gamma_X"]),   0.10)
})

test_that("simulate_counterfactual returns expected shapes", {
  sim <- simulate_tvhte(N = 300, T = 5, t0 = 3, J = 2,
                        beta = 0.4,
                        feedback_gamma = c(0.2, 0.3, 0.5, 0.4),
                        seed = 53)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X,
               compute_se = FALSE)
  fb  <- fit_feedback(sim$Y, sim$Y0, sim$X, sim$X0)

  cf <- simulate_counterfactual(fit, fb, t0_star = 4, N_star = 200, seed = 54)
  expect_equal(dim(cf$Y), c(200, 5))
  expect_equal(dim(cf$X), c(200, 5, 1))
  expect_equal(dim(cf$lambda), c(200, 2))
  expect_true(all(is.finite(cf$Y)))
})

test_that("counterfactual delays response when treatment timing is pushed back", {
  sim <- simulate_tvhte(N = 400, T = 6, t0 = 3, J = 2,
                        rho_Y = 0.3, rho_delta = 0.5, beta = 0.4,
                        feedback_gamma = c(0.2, 0.3, 0.5, 0.4),
                        seed = 55)
  fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X,
               compute_se = FALSE)
  fb  <- fit_feedback(sim$Y, sim$Y0, sim$X, sim$X0)
  cf_early <- simulate_counterfactual(fit, fb, t0_star = 3,
                                       N_star = 300, seed = 60)
  cf_late  <- simulate_counterfactual(fit, fb, t0_star = 5,
                                       N_star = 300, seed = 60)
  expect_gt(mean(cf_early$Y[, 4]) - mean(cf_late$Y[, 4]), 0.2)
})
