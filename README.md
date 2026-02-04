# bayessurvival <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Jamie-Wilson-UL/bayessurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jamie-Wilson-UL/bayessurvival/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The **bayessurvival** package provides Bayesian multiple imputation for right-censored survival data. It generates posterior draws for censored observations, allowing you to create completed datasets and apply standard survival summaries or downstream analyses with proper uncertainty accounting.

It supports:
- **Parametric** models (Weibull, exponential, lognormal) via Stan
- **Nonparametric** LDDP imputation
- **Grouped** analyses for comparisons across strata

## Installation

### Development version (GitHub)

```r
# install.packages("devtools")
devtools::install_github("Jamie-Wilson-UL/bayessurvival")
```

### Stan backend (rstan)

Parametric MCMC uses `rstan`. You will need a working C++ toolchain so Stan can compile the model on first use:

- **Windows**: install RTools
- **macOS**: `xcode-select --install`
- **Linux**: install `g++` and `make` (e.g. `build-essential`)

Recommended runtime settings:

```r
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## Quick Start

```r
library(bayessurvival)
library(survival)

fit <- impute(
  survival::lung,
  n_imputations = 5,
  model = "parametric",         # or "auto" / "nonparametric"
  distribution = "weibull"
)

print(fit)
plot(fit, type = "survival")

# Access completed datasets
first_dataset <- complete(fit, dataset = 1)
```

## One-Function Workflow

In most cases you can use `impute()` directly:

```r
fit <- impute(
  survival::lung,
  n_imputations = 5,
  distribution = "weibull",
  time_unit = "days"
)
```

`impute()` will try to detect your time and status columns automatically. If your column names are unusual, specify them explicitly:

```r
fit <- impute(
  survival::lung,
  time = "time",
  status = "status",
  n_imputations = 5,
  distribution = "lognormal"
)
```

Switching engines is just a change of arguments:

```r
fit_np  <- impute(survival::lung, model = "nonparametric", n_imputations = 5)
fit_grp <- impute(survival::lung, groups = "sex", distribution = "weibull", n_imputations = 5)
```

If you want direct access, the underlying functions are available:
- `bayesian_impute()` (parametric)
- `bayes_np_impute()` (nonparametric)

## Looking at Results

```r
print(fit)
```

Useful plots:

```r
plot(fit, type = "survival")                  # observed vs imputed survival curves
plot(fit, type = "completed_dataset_summary") # overview of one completed dataset
plot(fit, type = "boxplots_comparison")       # compare imputed datasets
plot(fit, type = "density")                   # logspline density comparison
```

Diagnostics (requires `bayesplot`):

```r
plot(fit, type = "posterior")
plot(fit, type = "trace")
plot(fit, type = "pairs")
```

Grouped fits use the same plotting API:

```r
print(fit_grp)
plot(fit_grp, type = "survival")
plot(fit_grp, type = "group_comparison")
```

## Working with Completed Datasets

```r
complete(fit, dataset = 1)      # a single completed dataset
all_sets <- complete(fit)       # list of datasets
long_fmt <- complete(fit, format = "long")
```

Completed datasets include transparency columns (`original_time`, `original_status`, `was_censored`).

For grouped analyses:

```r
combined <- complete(fit_grp, dataset = 1, groups = "combined")
separate <- complete(fit_grp, dataset = 1, groups = "separate")
```

## Exporting

```r
export(fit, "weibull_results", format = "csv")
export(fit, "weibull_results", format = "rds")
```

```r
export(fit, "clean_results", format = "csv", include_original = FALSE)
export(fit_grp, "groups", format = "csv", groups = "separate")
```

When exporting multiple datasets to RDS, set `combine = FALSE` to write one file per dataset.

## Advanced Usage

Custom priors and MCMC options are supported:

```r
my_priors <- get_default_priors("weibull")
my_priors$mu_log_shape <- 0.2
my_priors$sd_log_shape <- 0.4

custom_fit <- impute(
  survival::lung,
  time = "time",
  status = "status",
  distribution = "weibull",
  n_imputations = 20,
  priors = my_priors,
  mcmc_options = list(
    chains = 4,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.9,
    max_treedepth = 12
  )
)
```

Nonparametric runs accept `mcmc = list(nburn = ..., nsave = ..., nskip = ...)` in the same way.

## Troubleshooting

- **Columns not found**: specify `time` and `status` explicitly.
- **Stan compilation errors**: confirm C++ toolchain is installed (RTools / Xcode CLT / build-essential).
- **MCMC warnings**: increase iterations or chains.
- **Performance**: use fewer imputations or a simpler distribution (e.g., exponential).

## Documentation

- User guide: `docs/USER_GUIDE.md`
- Function reference: `?impute`, `?bayesian_impute`

## License

GPL (>= 2). See `LICENSE.md` for the full text.
