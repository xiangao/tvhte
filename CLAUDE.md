# tvhte — project notes for Claude

R package implementing Botosaru & Liu (2025, arXiv:2509.13698) and
(2026, *AEA P&P* 116: 70–74): a dynamic panel with correlated random
coefficients `(α_i, δ_{i0})`, AR(1) event-time effects `δ_{ij}`, and
optional homogeneous covariate feedback. Two-step QMLE +
Gaussian-conjugate empirical Bayes.

Companion Julia port: [`TVHTE.jl`](https://github.com/xiangao/TVHTE.jl).
Mirror changes both ways.

## What’s where

- `R/design.R` — `.unit_design()` builds the per-unit linear-Gaussian
  system matrices `(A, M, c, W, AMc, A1, AMW)`. `.unit_moments()` forms
  the marginal MVN mean and covariance of `Y_i | Y_{i0}` with the latent
  `(α_i, δ_{i0}, ε_i)` integrated out analytically (no Kalman filter —
  closed-form thanks to linear Gaussian).
- `R/simulate.R` —
  [`simulate_tvhte()`](https://xiangao.github.io/tvhte/reference/simulate_tvhte.md).
  Accepts staggered `t0` (length-N with `Inf` for never-treated),
  optional covariates via `beta`, optional feedback DGP via
  `feedback_gamma`.
- `R/tvhte.R` — main two-step estimator. Step 1 is BFGS on the marginal
  MVN log-likelihood; step 2 is closed-form Gaussian posterior (Tweedie
  collapses to standard conjugacy). Designs are cached per unique cohort
  for speed under staggered timing.
- `R/feedback.R` —
  [`fit_feedback()`](https://xiangao.github.io/tvhte/reference/fit_feedback.md)
  (pooled OLS for `X` dynamics) and
  [`simulate_counterfactual()`](https://xiangao.github.io/tvhte/reference/simulate_counterfactual.md)
  (Botosaru-Liu 2026 Algorithm 1).
- `inst/sims/mc_validation.R` — **parallel** MC harness using
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html). DGP
  1 N=500 results in `findings.md`; N=1500 scenarios are opt-in (slow
  even in parallel due to numerical Hessian step).
- `vignettes/illustrative.Rmd` — end-to-end story on a synthetic DGP
  with mixed cohorts, feedback, and direct/indirect counterfactual
  decomposition.

## Performance notes

- The numerical Hessian step in `optim(..., hessian = TRUE)` roughly
  doubles per-fit cost. Use `tvhte(..., compute_se = FALSE)` for tight
  loops (MC, vignette rendering).
- For staggered designs, build_ds() caches one `.unit_design` per unique
  cohort. This matters when many units share a cohort.

## Estimand pin / known issue from BL 2025

The QMLE under Gaussian working assumption on λ is consistent for common
parameters even when the true λ distribution is non-Gaussian. The
Hessian-based SE is NOT automatically robust — a sandwich SE is on the
roadmap (Phase 4 follow-up).

## Docs

pkgdown deploys to `gh-pages` on every push. Live at
<https://xiangao.github.io/tvhte/>. The function reference pages show
`@examples` blocks (one per exported function). Index page lists
functions only — examples live on per-function pages.

## Mirroring Julia changes

When the math changes: 1. Edit R first (faster iteration). 2. Mirror in
`TVHTE.jl/src/*.jl` keeping the same NamedTuple shape. 3. Verify both
pass their respective tests.
