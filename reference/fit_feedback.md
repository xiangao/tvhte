# Estimate the homogeneous covariate feedback process

Implements the covariate-feedback estimation step of Botosaru and Liu
(2026): under the homogeneous-feedback assumption that \$\$f_t(X\_{it}
\| I_i^{t-1}, \lambda_i) = f_t(X\_{it} \| I_i^{t-1})\$\$ (the covariate
adjustment rule is the same across units conditional on observable
history), the feedback process can be estimated separately from the
structural outcome model. The likelihood factors into a structural piece
(for `Y | X, history`) and a feedback piece (for `X | history`); this
function fits the latter.

## Usage

``` r
fit_feedback(Y, Y0, X, X0)
```

## Arguments

- Y:

  N x T outcome matrix (post-baseline).

- Y0:

  Length-N baseline outcome.

- X:

  N x T x 1 covariate array. (Multi-covariate feedback is a future
  extension; pass a single covariate for now.)

- X0:

  Length-N baseline covariate value, used as the lag at t = 1.

## Value

A list of class `"tvhte_feedback"` with:

- coef:

  Named numeric: `intercept, gamma_Y, gamma_X`.

- sigma_eta:

  Residual SD of the feedback equation.

- fit:

  The underlying `lm` object.

## Details

Phase 5 scope: a single covariate (`K = 1`), linear AR(1) feedback with
one lag of `Y` and one lag of `X`: \$\$X\_{it} = \gamma_0 + \gamma_Y
Y\_{i,t-1} + \gamma_X X\_{i,t-1} + \eta\_{it}.\$\$ Fit by pooled OLS
across unit-time observations with `t >= 2`.

## References

Botosaru, Irene and Laura Liu (2026). "Event Studies with Feedback." AEA
Papers and Proceedings 116: 70-74.
