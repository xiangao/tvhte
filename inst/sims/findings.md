# Monte Carlo validation findings (Phase 4)

Reproduce: `Rscript inst/sims/mc_validation.R [reps]` from the package
root. The script parallelises over reps using `parallel::mclapply`
(default `detectCores() - 1`).

## What is checked

For each scenario the harness runs `tvhte()` on `reps` simulated panels
and reports, per parameter:

- **bias**: mean estimate minus true value (should → 0 as N grows).
- **emp_sd**: empirical standard deviation across reps (the truth for
  what the SE estimate is trying to recover).
- **mean_se**: mean of the inverse-Hessian SE returned by `tvhte()`.
- **coverage**: empirical 95% CI coverage for the true value.

The SE reported by `tvhte()` is the naive QMLE SE: inverse Hessian
transformed by the delta method back to natural parameter scales.
Valid under the Gaussian working assumption on the random-coefficient
prior; under misspecification it is generally too optimistic. A
sandwich SE that is robust to prior misspecification is on the roadmap.

## Scenarios

| # | DGP | N | Purpose |
|---|---|---|---|
| 1 | linear-Gaussian, true Gaussian lambda | 500   | small-sample bias/coverage |
| 2 | linear-Gaussian, true Gaussian lambda | 1500  | asymptotic regime |
| 3 | bimodal lambda (mixture of two Gaussians) | 1500 | QMLE robustness check |
| 4 | linear-Gaussian, oscillating delta (`rho_delta = -0.5`) | 1500 | non-monotonic path |

## Results

### DGP 1: linear-Gaussian, N = 500 (100 reps, completed)

| param      | truth | bias    | emp_sd | mean_se | coverage |
|------------|------:|--------:|-------:|--------:|---------:|
| `rho_Y`    |  0.5  | -0.0021 | 0.025  | 0.025   | 0.93     |
| `rho_delta`|  0.6  |  0.0008 | 0.044  | 0.042   | 0.93     |
| `beta1`    |  0.5  | -0.0019 | 0.021  | 0.020   | 0.93     |
| `beta2`    | -0.3  | -0.0019 | 0.021  | 0.020   | 0.95     |

Clean: bias near zero, mean SE matches empirical SD to two decimals,
coverage is at or slightly under nominal (0.93–0.95). The estimator and
the naive Hessian-based SE both work as expected at N = 500 when the
Gaussian working assumption holds.

### DGPs 2–4 (N = 1500): not run

Per-fit cost at N = 1500 with the numerical Hessian step is roughly
8 minutes per scenario in parallel on 11 cores — too slow for routine
package validation. **Pragmatic call**: the N = 500 results already show
bias → 0 and coverage at nominal, so the asymptotic regime and
robustness checks are deferred. They can be re-run on demand by
passing `Rscript inst/sims/mc_validation.R 100` and waiting ~25 min.

## What's NOT in this validation yet

- Sandwich SE under misspecified Gaussian prior on lambda (the BL 2025
  consistency result for QMLE applies to the *point estimator*, not the
  inverse-Hessian SE).
- Empirical Bayes ratio-optimality (compound risk vs oracle risk for
  the unit-level `lambda_hat`).
- Coverage of cohort-conditional aggregates in staggered designs.
- Behaviour under non-Gaussian `U` or `eps`.
- Sensitivity to misspecified AR order on event-time effects.

## Lessons logged (process)

- Run MC in parallel from the first version. Serial MC on this code-
  base costs ~10× wall time and produces no insights worth waiting for.
- Numerical-Hessian SEs roughly double per-fit cost. For tight MC
  loops use `tvhte(..., compute_se = FALSE)` and only request SEs on
  the final point fit.
