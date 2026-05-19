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
  feedback_gamma = NULL,
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

- feedback_gamma:

  Optional list with components
  `c(intercept, gamma_Y, gamma_X, sigma_eta)` of length 4 enabling a
  single-covariate feedback DGP per Botosaru-Liu (2026): instead of iid
  X, generate \$\$X\_{it} = gamma_0 + gamma_Y Y\_{i,t-1} + gamma_X
  X\_{i,t-1} + \eta\_{it}.\$\$ Requires `length(beta) == 1`. Default
  `NULL` (no feedback).

- seed:

  Optional integer seed.

## Value

A list with `Y` (an N x T matrix), `Y0` (length-N baseline vector), `t0`
(scalar), `J`, `lambda` (N x 2 matrix of true alpha and delta_i0), and
`delta` (N x (J+1) matrix of true event-time effects).

## Details

All units share a common treatment period `t0` and event-time window
`0:J` (Phase 1 scope).

## Examples

``` r
# Common adoption timing
sim <- simulate_tvhte(N = 200, T = 5, t0 = 3, J = 2,
                      rho_Y = 0.4, rho_delta = 0.6, seed = 1)
dim(sim$Y); head(sim$Y0)
#> [1] 200   5
#> [1] -1.08690882 -1.82608301  0.99528181 -0.01186178 -0.59962839 -0.17794799

# Staggered adoption with Inf marking never-treated units
set.seed(2)
cohorts <- sample(c(3, 5, Inf), 200, replace = TRUE)
sim2 <- simulate_tvhte(N = 200, T = 7, t0 = cohorts, J = 2, seed = 2)
table(sim2$t0)
#> 
#>   3   5 Inf 
#>  71  69  60 

# With a covariate that evolves endogenously (Botosaru-Liu 2026 feedback)
sim3 <- simulate_tvhte(N = 200, T = 5, t0 = 3, J = 2,
                       beta = 0.4,
                       feedback_gamma = c(0.2, 0.3, 0.5, 0.4),
                       seed = 3)
head(sim3$X[, , 1])
#>             [,1]        [,2]         [,3]       [,4]       [,5]
#> [1,]  0.49404799 -0.40746550 -0.798322470  0.8795465  1.8417463
#> [2,]  0.07982317  0.54977545  0.681280804  0.6314203  1.3619035
#> [3,] -0.04297785  1.19607501  0.962776925  1.1901954  1.0792305
#> [4,] -0.19603226 -0.63674908 -0.071879100  0.5684659  1.0673672
#> [5,]  0.50722946  0.15350061  0.140027557 -0.2743175 -0.3239418
#> [6,] -0.50478521  0.07483239  0.004695086  0.8352969  2.7102571
```
