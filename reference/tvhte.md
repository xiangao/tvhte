# Two-step TV-HTE estimator (Botosaru and Liu 2025)

Fits the linear-Gaussian dynamic-panel model with correlated random
coefficients \$\$Y\_{it} = \rho_Y Y\_{i,t-1} + \alpha_i + \sum_j
D\_{it}^j \delta\_{ij} + U\_{it}\$\$ \$\$\delta\_{ij} = \rho\_\delta
\delta\_{i,j-1} + \varepsilon\_{ij},\\ j \geq 1\$\$ under common
treatment timing. The latent `lambda_i = (alpha_i, delta_{i0})` is
integrated out under a Gaussian working assumption; in step 1 we
estimate the common parameters (and the lambda-prior parameters) by QMLE
on the resulting marginal multivariate-normal log-likelihood; in step 2
we form the closed-form Gaussian posterior mean of
`(alpha_i, delta_{i0})` for each unit. The Botosaru-Liu (2025) paper
shows that the QMLE remains consistent under misspecification of the
Gaussian working assumption on lambda.

## Usage

``` r
tvhte(
  Y,
  Y0,
  t0,
  J,
  X = NULL,
  init = NULL,
  control = list(maxit = 500),
  compute_se = TRUE
)
```

## Arguments

- Y:

  N x T matrix of post-baseline outcomes.

- Y0:

  Length-N vector of baseline outcomes `Y_{i,0}`.

- t0:

  Treatment period. Either a scalar (common timing) or a
  length-`nrow(Y)` vector with per-unit cohorts (staggered). Use `Inf`
  for never-treated units.

- J:

  Maximum event time observed in-window (in `0:(T - t0)`).

- X:

  Optional `N x T x K` array of strictly exogenous covariates. If
  supplied, the model adds `X_{it}'beta` to the outcome equation and
  estimates `beta` jointly with the other parameters.

- init:

  Optional named list of starting values for the optimizer (overrides
  any of the defaults computed from the data).

- control:

  Passed to `optim`.

## Value

A list of class `"tvhte"` with:

- theta:

  Estimated common parameters `c(rho_Y, rho_delta, sigma_U, sigma_eps)`.

- prior:

  Estimated Gaussian prior on `lambda_i`.

- loglik:

  Maximised marginal log-likelihood.

- lambda_hat:

  N x 2 matrix of posterior means of `(alpha_i, delta_{i0})`.

- delta_path:

  N x (J+1) matrix of posterior-mean event-time trajectories. Column
  `j+1` is `E[delta_{i,j} | data]`.

- convergence:

  `optim`'s return code; 0 means converged.

## References

Botosaru, Irene and Laura Liu (2025). "Time-Varying Heterogeneous
Treatment Effects in Event Studies." arXiv:2509.13698.

## Examples

``` r
# Small panel, common adoption, no covariates
sim <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 3,
                      rho_Y = 0.5, rho_delta = 0.7, seed = 1)
fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J,
             compute_se = FALSE)
print(fit)
#> Time-Varying Heterogeneous Treatment Effects (Botosaru-Liu 2025)
#>   N = 300 units; T = 6 periods; t0 = 3; J = 3
#>   log-likelihood = -2808.536   convergence = 0
#> 
#> Common parameters (theta):
#>   rho_Y      = 0.5217
#>   rho_delta  = 0.7198
#>   sigma_U    = 1.033
#>   sigma_eps  = 0.349
#> Prior on lambda_i = (alpha, delta_0):
#>   mu = (0.02566, 0.9811)
#>   sd = (0.4622, 0.157),  cor = 0.9997
#> 
#> Mean posterior event-time effects (across units):
#>   delta_0 = 0.9811   delta_1 = 0.7084   delta_2 = 0.5116   delta_3 = 0.3640

# With a strictly-exogenous covariate
sim2 <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 3,
                       beta = 0.5, seed = 2)
fit2 <- tvhte(sim2$Y, sim2$Y0, t0 = sim2$t0, J = sim2$J,
              X = sim2$X, compute_se = FALSE)
fit2$beta
#> [1] 0.4770541

# Staggered adoption with a never-treated comparison group
set.seed(3)
cohorts <- sample(c(3, 5, Inf), 400, replace = TRUE)
sim3 <- simulate_tvhte(N = 400, T = 7, t0 = cohorts, J = 2, seed = 3)
fit3 <- tvhte(sim3$Y, sim3$Y0, t0 = sim3$t0, J = sim3$J,
              compute_se = FALSE)
fit3$cohort_counts
#> t0
#>   3   5 Inf 
#> 140 138 122 
```
