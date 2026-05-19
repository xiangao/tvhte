# Package index

## Main estimator

Two-step QMLE + Tweedie/Gaussian empirical-Bayes estimator (Botosaru-Liu
2025). Click through for runnable examples.

- [`tvhte()`](https://xiangao.github.io/tvhte/reference/tvhte.md) :
  Two-step TV-HTE estimator (Botosaru and Liu 2025)

## Covariate feedback (Botosaru-Liu 2026)

Homogeneous-feedback factorisation and counterfactual decomposition into
direct vs indirect components.

- [`fit_feedback()`](https://xiangao.github.io/tvhte/reference/fit_feedback.md)
  : Estimate the homogeneous covariate feedback process
- [`simulate_counterfactual()`](https://xiangao.github.io/tvhte/reference/simulate_counterfactual.md)
  : Simulate counterfactual paths under tvhte + feedback fits

## Simulation

Data generator covering vanilla, staggered (with never-treated), and
feedback DGPs. Click through for example calls.

- [`simulate_tvhte()`](https://xiangao.github.io/tvhte/reference/simulate_tvhte.md)
  : Simulate from the Botosaru-Liu (2025) TV-HTE model
