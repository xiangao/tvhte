# Simulate from the Botosaru-Liu (2025) TV-HTE model

Draws a dataset of N units observed at times 0, 1, ..., T according to
\$\$Y\_{it} = \rho_Y Y\_{i,t-1} + \alpha_i + \sum\_{j} D\_{it}^j
\delta\_{ij} + U\_{it}\$\$ \$\$\delta\_{ij} = \rho\_\delta
\delta\_{i,j-1} + \varepsilon\_{ij}, \quad j \geq 1\$\$ with
`(alpha_i, delta_{i0})` jointly Gaussian.

## Usage

``` r
simulate_tvhte(
  N = 500,
  T = 6,
  t0 = 3,
  J = 3,
  rho_Y = 0.5,
  rho_delta = 0.7,
  sigma_U = 1,
  sigma_eps = 0.3,
  mu_alpha = 0,
  mu_delta0 = 1,
  sigma_alpha = 0.5,
  sigma_delta0 = 0.5,
  cor_alpha_delta = 0,
  Y0_mean = 0,
  Y0_sd = 1,
  beta = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of units.

- T:

  Number of post-baseline periods (i.e. observations at t = 1, ..., T).

- t0:

  Treatment period. Either a scalar (common timing) in `1:T`, `Inf`
  (entire sample never treated), or a length-N integer/numeric vector
  with per-unit cohorts. Use `Inf` entries for never-treated units in
  the staggered case.

- J:

  Maximum event time observed in-window.

- rho_Y:

  Outcome AR(1) coefficient.

- rho_delta:

  Event-time AR(1) coefficient.

- sigma_U, sigma_eps:

  Standard deviations of `U_{it}` and `eps_{ij}`.

- mu_alpha, mu_delta0, sigma_alpha, sigma_delta0, cor_alpha_delta:

  Parameters of the Gaussian prior on
  `lambda_i = (alpha_i, delta_{i0})`.

- Y0_mean, Y0_sd:

  Distribution of the baseline outcome `Y_{i,0}`.

- beta:

  Optional numeric vector of length `K` of true coefficients on strictly
  exogenous covariates. If supplied, an `N x T x K` array of covariates
  is generated (standard normal by default) and added to the outcome
  equation as `X_{it}'beta`. Default `NULL` (no covariates).

- seed:

  Optional integer seed.

## Value

A list with `Y` (an N x T matrix), `Y0` (length-N baseline vector), `t0`
(scalar), `J`, `lambda` (N x 2 matrix of true alpha and delta_i0), and
`delta` (N x (J+1) matrix of true event-time effects).

## Details

All units share a common treatment period `t0` and event-time window
`0:J` (Phase 1 scope).
