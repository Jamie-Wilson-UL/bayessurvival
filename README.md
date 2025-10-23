# bayessurvival <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Jamie-Wilson-UL/bayessurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jamie-Wilson-UL/bayessurvival/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The **bayessurvival** package provides Bayesian multiple imputation for right-censored survival data. It generates posterior draws for censored observations, allowing you to create completed datasets and apply standard survival summaries or downstream analyses with proper uncertainty accounting.

## Installation

### Development version

```r
# install.packages("devtools")
devtools::install_github("Jamie-Wilson-UL/bayessurvival")
```

### CmdStan requirement

Bayesian estimation is carried out via `cmdstanr`. If CmdStan is not available, you will be prompted to install it:

```r
cmdstanr::install_cmdstan()
```

## Quick Start

```r
library(bayessurvival)
library(survival)

# Prepare and impute in one step (automatic column detection)
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

For nonparametric LDDP imputation, set `model = "nonparametric"` (or use `distribution = "nonparametric"`).

## Documentation

- User guide and workflow notes: `USER_GUIDE.md`
- Function reference: `pkgdown` site (if available) or `?bayesian_impute`

## License

GPL (>= 2). See [LICENSE.md](LICENSE.md) for the full text. 