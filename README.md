<!-- README.md is generated from README.Rmd. Please edit that file -->



# SelectBoost.beta <img src="man/figures/logo.svg" align="right" width="200"/>

<!-- badges: start -->
<!-- [![DOI](https://img.shields.io/badge/doi-10.32614/CRAN.package.SelectBoost.beta-blue.svg)](https://doi.org/10.32614/CRAN.package.SelectBoost.beta) -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/SelectBoost.beta)](https://cran.r-project.org/package=SelectBoost.beta) -->
[![R-CMD-check](https://github.com/fbertran/SelectBoost.beta/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/SelectBoost.beta/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/fbertran/SelectBoost.beta/actions/workflows/rhub.yaml/badge.svg)](https://github.com/fbertran/SelectBoost.beta/actions/workflows/rhub.yaml)
<!-- badges: end -->


With the growth of big data, variable selection has become one of the major challenges in statistics. Although many methods have been proposed in the literature their performance in terms of recall and precision are limited in a context where the number of variables by far exceeds the number of observations or in a high correlated setting. 

Results: 

`SelectBoost.beta` brings the correlation-aware resampling strategy of the original SelectBoost package to **beta** regression by implementing an extension of the **SelectBoost** algorithm, F. Bertrand, I. Aouadi, N. Jung, R. Carapito, L. Vallat, S. Bahram, M. Maumy-Bertrand (2015) <https://doi.org/10.1093/bioinformatics/btaa855> and <https://doi.org/10.32614/CRAN.package.SelectBoost>.


It ships with:

- wrappers such as `betareg_step_aic()` and `betareg_glmnet()` that act as base
  selectors for beta-distributed outcomes, now including optional precision
  (phi) submodel search and observation weights;
- helper functions (`sb_normalize()`, `sb_group_variables()`,
  `sb_resample_groups()`, …) mirroring the core stages of SelectBoost; and
- the high-level `sb_beta()` driver that orchestrates normalisation,
  correlation analysis, grouped resampling and stability tallying in a single
  call.

Each resampling call returns per-group diagnostics (cached draws, observed
correlation summaries) and `sb_beta()` threads the same correlated surrogates
across all thresholds so cross-level comparisons remain aligned. Interval
responses are supported through the `interval` argument, which reuses the
`fastboost_interval()` logic directly inside `sb_beta()`.

The package is designed so that each stage of the workflow remains reusable on
its own. Users can plug in custom grouping strategies or selectors while still
benefiting from correlated resampling.



## Conference presentations

The SelectBoost4Beta approach was presented by Frédéric Bertrand and Myriam
Maumy at the Joint Statistical Meetings 2023 in Toronto ("Improving variable
selection in Beta regression models using correlated resampling") and at
BioC2023 in Boston ("SelectBoost4Beta: Improving variable selection in Beta
regression models"). Both communications highlighted how correlated resampling boosts
variable selection for Beta regression in high-dimensional, strongly correlated
settings.

## Installation

SelectBoost.beta is preparing for its first CRAN submission. Until it becomes
available there, install the development version from GitHub:


``` r
devtools::install_github("fbertran/SelectBoost.beta")
```

Once the package lands on CRAN, the usual `install.packages("SelectBoost.beta")`
command will work as expected.

The selectors rely on the `betareg`, `glmnet`, and `gamlss` ecosystems. These
packages will be pulled in automatically when installing from source.

## Quick start

Simulate a correlated design, run the manual SelectBoost steps with
`betareg_step_aic()`, and compute selection frequencies:


``` r
library(SelectBoost.beta)
set.seed(42)

sim <- simulation_DATA.beta(n = 150, p = 6, s = 3, beta_size = c(1, -0.8, 0.6))
X_norm <- sb_normalize(sim$X)
corr_mat <- sb_compute_corr(X_norm)
groups <- sb_group_variables(corr_mat, c0 = 0.6)
resamples <- sb_resample_groups(X_norm, groups, B = 50)
#> Error in group_draws[[i]]: subscript out of bounds
coef_path <- sb_apply_selector_manual(X_norm, resamples, sim$Y, betareg_step_aic)
#> Error in .betareg_step_engine(X = X, Y = Y, criterion_fun = function(fit) .aic_betareg(fit, : Initial betareg fit failed: Error in model.frame.default(terms(formula, lhs = lhs, rhs = rhs, data = data,  : 
#>   invalid type (closure) for variable '(weights)'
sel_freq <- sb_selection_frequency(coef_path, version = "glmnet")
sel_freq
#> x1 x2 x3 x4 x5 x6 
#>  1  1  1  0  0  0

attr(resamples, "diagnostics")
#> NULL
```

The `sb_beta()` wrapper performs the entire loop internally and returns a matrix
indexed by the correlation thresholds used during resampling:


``` r
sb <- sb_beta(sim$X, sim$Y, B = 50, step.num = 0.25,use.parallel = FALSE)
print(sb)
#> SelectBoost beta selection frequencies
#> Selector: betareg_step_aic
#> Resamples per threshold: 50
#> c0 grid: 1.000, 0.089, 0.059, 0.030, 0.000
#> Inner thresholds: 0.089, 0.059, 0.030
#>              x1   x2   x3   x4   x5   x6
#> c0 = 1.000 1.00 1.00 1.00 0.00 0.00 0.00
#> c0 = 0.089 0.22 0.18 0.10 0.18 0.18 0.16
#> c0 = 0.059 0.10 0.12 0.20 0.20 0.16 0.28
#> c0 = 0.030 0.16 0.20 0.12 0.08 0.16 0.14
#> c0 = 0.000 0.28 0.22 0.12 0.08 0.18 0.18
#> attr(,"c0.seq")
#> [1] 1.00000000 0.08894615 0.05949716 0.03010630 0.00000000
#> attr(,"steps.seq")
#> [1] 0.08894615 0.05949716 0.03010630
#> attr(,"B")
#> [1] 50
#> attr(,"selector")
#> [1] "betareg_step_aic"
```

The result stores the selector used, the number of resamples, and the
correlation thresholds in its attributes. Dedicated methods make these easier to
inspect programmatically:


``` r
summary(sb)
#> SelectBoost beta summary
#> Selector: betareg_step_aic
#> Resamples per threshold: 50
#> c0 grid: 1.000, 0.089, 0.059, 0.030, 0.000
#> Inner thresholds: 0.089, 0.059, 0.030
#> Top rows:
#>        c0 variable frequency
#> 1  1.0000       x1      1.00
#> 2  1.0000       x2      0.22
#> 3  1.0000       x3      0.10
#> 4  1.0000       x4      0.16
#> 5  1.0000       x5      0.28
#> 6  1.0000       x6      1.00
#> 7  0.0889       x1      0.18
#> 8  0.0889       x2      0.12
#> 9  0.0889       x3      0.20
#> 10 0.0889       x4      0.22
if (requireNamespace("ggplot2", quietly = TRUE)) {
  autoplot(sb)
}
#> Error in autoplot(sb): could not find function "autoplot"
```


``` r
attr(sb, "selector")
#> [1] "betareg_step_aic"
attr(sb, "c0.seq")
#> [1] 1.00000000 0.08894615 0.05949716 0.03010630 0.00000000
attr(sb, "resample_diagnostics")[[1]]
#> NULL
```


``` r
single <- compare_selectors_single(sim$X, sim$Y, include_enet = TRUE)
```


``` r
freq <- suppressWarnings(compare_selectors_bootstrap(
  sim$X, sim$Y, B = 100, include_enet = TRUE, seed = 321
))
head(freq)
#>    selector variable freq
#> x1      AIC       x1 1.00
#> x2      AIC       x2 1.00
#> x3      AIC       x3 1.00
#> x4      AIC       x4 0.27
#> x5      AIC       x5 0.14
#> x6      AIC       x6 0.19
```



``` r
plot_compare_coeff(single$table)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-8-1.png" alt="plot of chunk unnamed-chunk-8" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-8</p>
</div>



``` r
plot_compare_freq(freq)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-9-1.png" alt="plot of chunk unnamed-chunk-9" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-9</p>
</div>

### Interval outcomes

`sb_beta()` can draw pseudo-responses from observed intervals by supplying
`Y_low`, `Y_high`, and an `interval` mode:


``` r
interval_fit <- sb_beta(
  sim$X,
  Y_low = pmax(sim$Y - 0.05, 0),
  Y_high = pmin(sim$Y + 0.05, 1),
  interval = "uniform",
  B = 30,
  step.num = 0.5
)
#> Error in sb_beta(sim$X, Y_low = pmax(sim$Y - 0.05, 0), Y_high = pmin(sim$Y + : object 'y' not found
attr(interval_fit, "interval")
#> Error: object 'interval_fit' not found
head(attr(interval_fit, "resample_diagnostics")[[length(interval_fit)]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'head': object 'interval_fit' not found
```

### Parallel resampling

Setting `use.parallel = TRUE` instructs `sb_beta()` and `sb_resample_groups()` to
dispatch resamples and selector fits through
[`future.apply`](https://future.apply.futureverse.org/). Bring your own
`future::plan()` to select the desired backend (e.g. `multisession` on desktops):


``` r
future::plan(future::multisession)
sb_parallel <- sb_beta(sim$X, sim$Y, B = 50, step.num = 0.25, use.parallel = TRUE)
future::plan(future::sequential)
```

Refer to the vignettes for a more detailed walk-through of the workflow and the
pseudo-code underpinning the algorithms.
