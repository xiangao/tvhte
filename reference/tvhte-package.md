# tvhte: Time-Varying Heterogeneous Treatment Effects in Event Studies

Two-step semiparametric estimator for event studies with persistent
outcomes and unit-level heterogeneous, dynamic treatment effects,
following Botosaru and Liu (2025). Imposes an autoregressive structure
on event-time effects to keep dimensionality manageable. First step:
quasi-maximum-likelihood for common parameters under a Gaussian working
assumption on the random coefficients. Second step: empirical Bayes
(Tweedie / Gaussian conjugacy) for unit-specific treatment trajectories.
A separate option implements the homogeneous feedback extension of
Botosaru and Liu (2026), which decomposes dynamic responses into direct
and indirect (covariate-adjustment) components.

## See also

Useful links:

- <https://xiangao.github.io/tvhte/>

- <https://github.com/xiangao/tvhte>

- Report bugs at <https://github.com/xiangao/tvhte/issues>

## Author

**Maintainer**: Xiang Ao <xiangao@gmail.com>

Authors:

- Xiang Ao <xiangao@gmail.com>
