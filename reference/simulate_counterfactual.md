# Simulate counterfactual paths under tvhte + feedback fits

Implements Algorithm 1 of Botosaru and Liu (2026): given fits of the
outcome model (`tvhte`) and feedback model (`fit_feedback`), simulate
counterfactual joint trajectories `{Y_i^{T,*}, X_i^{T,*}}` under
alternative treatment timing or initial conditions.

## Usage

``` r
simulate_counterfactual(
  fit,
  fit_fb,
  t0_star,
  N_star = 500,
  Y0_star = NULL,
  X0_star = NULL,
  seed = NULL
)
```

## Arguments

- fit:

  A `tvhte` object.

- fit_fb:

  A `tvhte_feedback` object.

- t0_star:

  Counterfactual treatment timing (scalar applied to all simulated
  units, or length-`N_star` vector). Use `Inf` for never- treated
  counterfactuals.

- N_star:

  Number of counterfactual units to simulate.

- Y0_star, X0_star:

  Counterfactual baseline values (length-`N_star` each). Default draws
  standard normals.

- seed:

  Optional RNG seed.

## Value

A list with `Y` (`N_star x T` matrix of counterfactual outcomes), `X`
(`N_star x T x 1` array), and `lambda` (`N_star x 2` matrix of drawn
latent effects).

## Details

Latent `lambda_i = (alpha_i, delta_{i0})` is drawn from the estimated
Gaussian prior; event-time effects evolve via the AR(1) on `delta_{ij}`;
covariates evolve via the estimated feedback process; outcomes
accumulate the direct effect (`delta_{ij}`) plus the indirect effect
(`X_it' beta`).

## Examples

``` r
# Fit the model, fit the feedback, then simulate counterfactual paths
# under shifted treatment timing.
sim <- simulate_tvhte(N = 300, T = 6, t0 = 3, J = 2, beta = 0.4,
                      feedback_gamma = c(0.2, 0.3, 0.5, 0.4), seed = 1)
fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J, X = sim$X,
             compute_se = FALSE)
fb  <- fit_feedback(sim$Y, sim$Y0, sim$X, sim$X0)

# Counterfactual: treatment delayed from t=3 to t=5
cf <- simulate_counterfactual(fit, fb, t0_star = 5,
                               N_star = 200, seed = 99)
colMeans(cf$Y)   # average counterfactual outcome path
#> [1] 0.1401096 0.2805138 0.3674985 0.3779429 1.6185435 1.9285446
```
