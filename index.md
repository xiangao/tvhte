# tvhte

[![pkgdown](https://img.shields.io/badge/pkgdown-site-blue.svg)](https://xiangao.github.io/tvhte/)

`tvhte` estimates time-varying heterogeneous treatment effects in event
studies, following Botosaru and Liu (2025, 2026).

## What this package does

Standard event-study regressions with unit and time fixed effects leave
little room for residual serial dependence. When outcomes are
persistent, the event-time dummies can pick up persistence as well as
the treatment effect. That can look like pre-trends or biased
post-treatment effects.

`tvhte` implements Botosaru & Liu’s two-step semiparametric estimator:

1.  **Marginal-likelihood QMLE** for common parameters
    `(ρ_Y, ρ_δ, σ_U, σ_ε)` and the latent prior, integrating out
    unit-level random coefficients `(α_i, δ_{i0})` under a Gaussian
    working assumption. Consistent under misspecification of the prior.
2.  **Tweedie / Gaussian-conjugate empirical Bayes** for posterior means
    of the unit-specific treatment trajectories `{δ_{i,j}}`, achieving
    asymptotic ratio optimality.

The dimensionality reduction comes from imposing an AR(1) on event-time
effects, `δ_{ij} = ρ_δ δ_{i,j-1} + ε_{ij}`, with `δ_{i0}` free.

## Status

The current implementation is still staged. The linear-Gaussian
common-adoption case is the most mature part. The planned pieces are:

- Phase 2: strictly exogenous covariates `X_{it}'β`.
- Phase 3: staggered adoption with cohort-specific conditioning.
- Phase 4: full MC validation suite.
- Phase 5: feedback / endogenous covariate adjustment (Botosaru-Liu
  2026).
- Phase 6: Great Recession county-unemployment replication vignette.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("xiangao/tvhte")
```

## Minimal example

``` r

library(tvhte)

# Simulate from the model
sim <- simulate_tvhte(N = 800, T = 6, t0 = 3, J = 3,
                      rho_Y = 0.5, rho_delta = 0.7,
                      sigma_alpha = 0.5, sigma_delta0 = 0.4,
                      cor_alpha_delta = 0.3, seed = 1)

# Fit
fit <- tvhte(sim$Y, sim$Y0, t0 = sim$t0, J = sim$J)
print(fit)
# rho_Y      = 0.5112
# rho_delta  = 0.7122
# sigma_U    = 1.007
# ...
```

## Documentation & vignettes

Full documentation: **<https://xiangao.github.io/tvhte/>**

| Page | Description |
|----|----|
| [Home](https://xiangao.github.io/tvhte/) | Overview, install, motivation |
| [Illustrative vignette](https://xiangao.github.io/tvhte/articles/illustrative.html) | End-to-end walk-through: simulate, fit, posterior trajectories, counterfactual direct/indirect decomposition, with plots |
| [`tvhte()`](https://xiangao.github.io/tvhte/reference/tvhte.html) | Main estimator (QMLE + Tweedie EB). Runnable examples: common adoption, with covariate, staggered with never-treated |
| [`simulate_tvhte()`](https://xiangao.github.io/tvhte/reference/simulate_tvhte.html) | Data generator (vanilla DGP and feedback DGP) |
| [`fit_feedback()`](https://xiangao.github.io/tvhte/reference/fit_feedback.html) | Botosaru-Liu 2026 homogeneous-feedback estimator |
| [`simulate_counterfactual()`](https://xiangao.github.io/tvhte/reference/simulate_counterfactual.html) | Algorithm 1 of BL 2026: joint counterfactual `(Y*, X*)` |
| [Reference index](https://xiangao.github.io/tvhte/reference/index.html) | All functions on one page |

A Julia port is available at
[TVHTE.jl](https://github.com/xiangao/TVHTE.jl)
([docs](https://xiangao.github.io/TVHTE.jl/)).

## References

- Botosaru, Irene and Laura Liu (2025). “Time-Varying Heterogeneous
  Treatment Effects in Event Studies.”
  [arXiv:2509.13698](https://arxiv.org/abs/2509.13698).
- Botosaru, Irene and Laura Liu (2026). “Event Studies with Feedback.”
  *AEA Papers and Proceedings* 116: 70–74.

## License

MIT
