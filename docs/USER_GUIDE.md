# bayessurvival User Guide

---

## Getting Started

Install from GitHub and load the library:

```r
install.packages("devtools")
devtools::install_github("Jamie-Wilson-UL/bayessurvival")

library(bayessurvival)
```

If you plan to use the parametric engines (Weibull, exponential, lognormal), you'll also need CmdStan the first time:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

**Important**: Installing CmdStan requires a C++ toolchain. If you encounter errors during `cmdstanr::install_cmdstan()`, you may need to install additional development tools:

- **Windows**: Install [RTools](https://cran.r-project.org/bin/windows/Rtools/) (includes C++ compiler)
- **macOS**: Install Xcode Command Line Tools: `xcode-select --install` in Terminal
- **Linux**: Install `g++` and `make` (usually pre-installed): `sudo apt-get install build-essential` (Ubuntu/Debian)

You can check if your toolchain is ready with:
```r
cmdstanr::check_cmdstan_toolchain()
```

The nonparametric LDDP engine doesn't rely on CmdStan, so it will work even if you skip that step.

Throughout the examples below, we use the classic `survival::lung` dataset:

```r
library(survival)
lung <- survival::lung
```

---

## One-Function Workflow

For most situations, you can just use `impute()`:

```r
fit <- impute(
  lung,
  n_imputations = 5,
  distribution = "weibull",    # choose "weibull", "exponential", or "lognormal"
  time_unit = "days"
)
```

By default, `impute()` tries to detect the time and status columns (looking for common names like `time`, `status`, `duration`, etc.). If your column names are unusual it’s advisable to specify them explicitly:

```r
fit <- impute(
  lung,
  time = "time",           # you can also pass indices, e.g. time = 3
  status = "status",
                                 n_imputations = 5,
  distribution = "lognormal"
)
```

Switching engines (distributions/groups) is just a change of arguments:

```r
fit_np  <- impute(lung, model = "nonparametric", n_imputations = 5)
fit_grp <- impute(lung, groups = "sex", distribution = "weibull", n_imputations = 5)
```

If you need direct access to the underlying functions, they are available as `bayesian_impute()` (parametric) and `bayes_np_impute()` (nonparametric). They behave exactly the same as `impute()` once you pass columns or `prepare_survival_data(...)` output.

---

## Looking at the Results

The print method gives a basic statistical summary with sensible next steps:

```r
print(fit)
```

You can sanity-check your imputations with the plotting helpers:

```r
plot(fit, type = "survival")                  # observed vs imputed survival curves
plot(fit, type = "completed_dataset_summary") # overview of a single (random) completed dataset
plot(fit, type = "boxplots_comparison")       # compare all imputed datasets side-by-side
plot(fit, type = "density")                   # logspline density comparison
```

When you want to dive deeper into model diagnostics, these are available too:

```r
plot(fit, type = "posterior")  # (random) posterior distributions
plot(fit, type = "trace")      # MCMC convergence (requires bayesplot package)
plot(fit, type = "pairs")      # parameter relationships (requires bayesplot package)
```

Note: The `trace` and `pairs` plot types rely on the `bayesplot` package. If you see a message asking you to install it, just run `install.packages(c("bayesplot", "posterior"))` and the plots will work.

All of the group-specific objects (`impute(..., groups = ...)`) respond to the same plotting API. For example:

```r
print(fit_grp)
plot(fit_grp, type = "survival")
plot(fit_grp, type = "group_comparison")
plot(fit_grp, type = "boxplots_comparison")
```

---

## Working with Completed Datasets

The package stores completed datasets in a way that’s deliberately transparent. You can pull out one dataset or many, in either wide, long or list form:

```r
complete(fit, dataset = 1)            # a single completed dataset
all_sets <- complete(fit)             # list of datasets
long_fmt <- complete(fit, format = "long")
```

Each dataset carries the original time/status columns alongside the imputed versions (`original_time`, `original_status`, `was_censored`), making it easy to see exactly what changed.

For grouped analyses you can choose to combine groups or keep them separate:

```r
combined <- complete(fit_grp, dataset = 1, groups = "combined")
separate <- complete(fit_grp, dataset = 1, groups = "separate")
```

---

## Exporting

When you’re satisfied with a fit, export the data in the format that suits your downstream work:

```r
export(fit, "weibull_results", format = "csv")
export(fit, "weibull_results", format = "rds")
```

There are helper arguments if you want to drop the transparency columns or focus on particular datasets:

```r
export(fit, "clean_results", format = "csv", include_original = FALSE)
export(fit_grp, "groups", format = "csv", groups = "separate")
```

At the moment only CSV and RDS exports are enabled in the package. (Support for Excel/Stata/SPSS is planned, but not yet implemented.)

---

## Advanced Usage

The defaults should be sufficient for most applications, but it is also possible to specify custom priors and MCMC settings. A good pattern is to start from the defaults and tweak the pieces you care about:

```r
# Start from the default Weibull priors and modify them
my_priors <- get_default_priors("weibull")
my_priors$mu_log_shape <- 0.2
my_priors$sd_log_shape <- 0.4
my_priors$mu_log_scale <- 0.5
my_priors$sd_log_scale <- 0.3

custom_fit <- impute(
  lung,
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
  ),
  verbose = TRUE
)
```

You can do the same for other distributions (`get_default_priors("lognormal")`, etc.). Nonparametric runs accept `mcmc = list(nburn = ..., nsave = ..., nskip = ...)` in the same way.

---

## Troubleshooting

A few quick checks:

- **Columns not found**: inspect with `colnames(lung)` and specify `time = "..."`, `status = "..."` explicitly (indices work too).
- **CmdStan not installed**: run `cmdstanr::install_cmdstan()` or set `model = "nonparametric"` temporarily.
- **MCMC diagnostics**: call `plot(fit, type = "trace")` to inspect convergence; increase iterations if needed.
- **Performance**: use fewer imputations (e.g., 5 instead of 50) or choose a faster distribution (`"exponential"`).
- **Reusing draws**: if you want 100 datasets but already ran 20, use `generate_complete_datasets()` on the stored posterior draws rather than refitting.

A quick sanity check to confirm your data meet the minimum requirements:

```r
stopifnot(
  all(lung$time > 0),
  all(lung$status %in% c(0, 1)),
  !any(is.na(lung$time)),
  !any(is.na(lung$status))
)
```

---

## Summary

- `impute()` is the default entry point.
- Parametric vs nonparametric vs grouped behaviour is controlled with `distribution`, `model`, and `groups` arguments.
- `complete()` and `export()` once you have a fit.
- Plots provide quality checks and diagnostics.
- Custom priors/MCMC settings are available if desired.