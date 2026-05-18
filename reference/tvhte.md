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
tvhte(Y, Y0, t0, J, X = NULL, init = NULL, control = list(maxit = 500))
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
