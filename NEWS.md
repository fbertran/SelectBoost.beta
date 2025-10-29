# SelectBoost.beta 0.4.5
- Added a pseudo-code vignette, refreshed the README with workflow details, and
  expanded the unit test suite for the new helpers.
- Fixed code and descriptions to get rid of notes during CRAN checks.
- Enabled optional parallel resampling in `sb_beta()`/`sb_resample_groups()` via
  `future.apply`, added S3 print/summary/autoplot helpers for `sb_beta`
  results, and documented the new behaviour in the README.

- Extended the stepwise beta selectors to handle observation weights and
  precision-submodel search, exposing precision coefficients in the returned
  paths.
- Added reproducible resampling caches and quality diagnostics to
  `sb_resample_groups()`/`sb_beta()`, including interval-response support that
  reuses pseudo-responses across correlation thresholds.
- Documented interval workflows more prominently by adding `sb_beta_interval()`,
  expanding the README/CRAN vignette guidance for selector choice and interval
  stability, and clarifying comparison-helper outputs and response squeezing.
  
# SelectBoost.beta 0.4.4
- Add **NEWS.md** documenting development history.

# SelectBoost.beta 0.4.3
- Added `sb_beta()` to run the full SelectBoost correlated-resampling loop with
  beta-regression selectors, plus a vignette illustrating the workflow.
- Added vignette section demonstrating the extended simulator and interval selection.

# SelectBoost.beta 0.4.2
- New **`simulation_DATA()`** to generate **interval-valued** Beta-regression data:
  - `interval = "jitter"` (symmetric) or `"quantile"` (Beta quantile intervals).
  - Works with `fastboost_interval()`; added a small vignette and unit test.
  - Supports *mixed mechanisms* (row-wise jitter vs quantile), *asymmetric* jitter widths (`delta_low`/`delta_high`),
    *asymmetric* quantile coverage (`alpha_low`/`alpha_high`), covariate-driven parameters (accept functions of `(mu, X)`),
    and optional **missing bounds** per row (`na_rate`, `na_side`).

# SelectBoost.beta 0.4.1
- **Comparison helpers** and **visualizations**:
  - `compare_selectors_single()`, `compare_selectors_bootstrap()` to run all selectors (AIC/BIC/AICc, GAMLSS LASSO/ENet*, GLMNET) and compute **selection frequencies**.
  - `plot_compare_coeff()`, `plot_compare_freq()` heatmaps to compare selectors **side by side**.
  - Vignette expanded to include simulated data and two real-ish datasets after scaling to (0,1).
- Kept `fastboost_interval()` (interval response stability selection), C++ IRLS speedups, and `prestandardize` option for `betareg_glmnet()`.
  - *ENet requires `gamlss.lasso` if installed.*

# SelectBoost.beta 0.4.0
- Initial beta-regression integration for SelectBoost:
  - Stepwise **AIC**, **BIC**, and **AICc** wrappers around `betareg`.
  - **GAMLSS** LASSO (`gamlss::ri`) and optional Elastic-Net (`gamlss.lasso::gnet`).
  - **Pure glmnet** (IRLS + weighted Gaussian ENet) with optional **prestandardize**.
  - `fastboost_interval()` prototype for interval responses.
  - First vignette + roxygenized docs.
