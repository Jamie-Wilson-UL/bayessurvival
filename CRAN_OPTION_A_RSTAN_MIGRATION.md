# Option A (CRAN): migrate Stan backend from `cmdstanr` to `rstan` / `rstantools`

This document maps the work required to implement “Option A”: keep the Stan-based parametric MCMC, but move from `cmdstanr` (CmdStan) to a CRAN-native Stan toolchain (`rstan`, ideally via `rstantools`) so CRAN can build/install the package from source without extra repositories or internet downloads.

## Why this is needed (current blockers)

In the current repo:

- `DESCRIPTION` declares `Additional_repositories: https://mc-stan.org/r-packages/` and `Suggests: cmdstanr` (`DESCRIPTION:45`, `DESCRIPTION:60`). CRAN generally requires that all dependencies are resolvable from CRAN and does not accept “extra repositories” metadata.
- The parametric engine compiles and runs Stan via `cmdstanr` (`R/stan-models.R`, `R/bayesian-impute.R`). This implies an external CmdStan toolchain that CRAN checks will not install/download.

Option A replaces this with `rstan` (on CRAN), optionally using `rstantools` so the Stan model compilation happens as part of `R CMD INSTALL` rather than at runtime.

## Current architecture (what depends on `cmdstanr`)

### Stan model compilation / caching
- `R/stan-models.R`
  - `compile_stan_models()` calls `cmdstanr::cmdstan_model()`
  - `get_stan_model()` caches compiled `CmdStanModel` objects in `.bayessurvival_env`
  - `stan_models_ready()` checks the cache

### Stan sampling and extraction
- `R/bayesian-impute.R`
  - `bayesian_impute_single()` calls `stan_model <- get_stan_model(distribution)` then `stan_fit <- stan_model$sample(...)`
  - Draws extraction is cmdstanr-shaped:
    - `posterior_samples <- stan_fit$draws(param_names)`
    - `extract_posterior_distributions()` uses `stan_fit$draws("t_imputed", format = "array")`
    - `compute_diagnostics()` uses `stan_fit$diagnostic_summary()` and `stan_fit$summary(variables = ...)`
    - `compute_waic_from_stan()` uses `stan_fit$draws()` and `stan_fit$draws(variables = ...)`
- `R/utils.R`
  - `extract_mcmc_summary()` uses `stan_fit$summary(...)` and `stan_fit$diagnostic_summary()`
  - `check_cmdstan()` and `check_required_packages()` refer to `cmdstanr`

### Tests / docs that mention CmdStan
- `README.md` instructs users to install CmdStan via `cmdstanr::install_cmdstan()`
- Multiple tests skip on `cmdstanr` (e.g., `tests/testthat/test-export.R`, `tests/testthat/test-complete.R`)

## Target architecture (Option A)

There are two realistic implementations under “Option A”:

### A1. `rstan` runtime compilation (smaller change, slower first-run)
Compile `.stan` to a `stanmodel` on first use via `rstan::stan_model()` (cache it), then sample via `rstan::sampling()`.

Pros:
- Minimal package build complexity (no generated C++ files committed).
- Lower upfront migration effort.

Cons:
- First-run compile can be slow for users.
- Still requires a working C++ toolchain at runtime.
- Not as “install-time compiled” as typical CRAN Stan packages.

### A2. `rstantools` install-time compilation (best CRAN story, more build work)
Use `rstantools` scaffolding so `.stan` files are translated/compiled as part of `R CMD INSTALL` (the common CRAN pattern for Stan packages).

Pros:
- Users do not compile Stan models at runtime (faster first use).
- Most aligned with “Stan-on-CRAN” conventions.

Cons:
- More invasive build-system changes and careful integration with existing `src/` Fortran/C.
- Higher risk of platform-specific build issues (Windows/macOS toolchains).

Recommendation: implement A1 first to validate correctness and refactor the API around an “rstan fit”, then move to A2 when stable.

## Work breakdown (high level)

### Phase 0 — CRAN hygiene prerequisites (do this regardless)
These issues are independent of the Stan backend but will block CRAN:
- Remove committed build artifacts:
  - `src/*.o`, `src/*.so` (these caused `R CMD check` warnings previously)
  - any `inst/stan/<executable>` (this repo currently has `inst/stan/weibull_imputation`)
- Remove stray files that could get bundled in the tarball:
  - `tests/.DS_Store`
  - root `lung_imputed_dataset_1.csv` (should not ship as package content unless intentional and documented)
- Remove `Additional_repositories` from `DESCRIPTION` and avoid non-CRAN dependencies.

### Phase 1 — Dependency + metadata changes (Stan backend swap)
Update package metadata to replace `cmdstanr` with `rstan`:
- `DESCRIPTION`
  - Remove `Additional_repositories`
  - Remove `cmdstanr` from `Suggests`
  - Add `rstan` (and potentially `StanHeaders`) as `Imports` (if parametric engine is core)
  - If implementing A2 (`rstantools`): add `LinkingTo`/`Imports` required by the scaffold (commonly `Rcpp`, `RcppEigen`, `StanHeaders`, `rstan`)
  - Update `SystemRequirements` to reflect a C++ toolchain requirement (no CmdStan download)
- `README.md`
  - Replace the CmdStan installation section with `rstan` toolchain guidance
- Man pages
  - Remove/replace `check_cmdstan` docs (either delete the function or replace with an `rstan`-focused check)

### Phase 2 — Replace model compilation/caching (`R/stan-models.R`)
Replace the cmdstanr-centric lifecycle with rstan equivalents:

New responsibilities for `get_stan_model()`:
- Return a cached `rstan::stanmodel`
- Compile on demand (A1) or return a prebuilt model object (A2)

Implementation notes:
- **A1:** `rstan::stan_model(file = system.file("stan", "<model>.stan", package = "bayessurvival"))`
- **A2:** `rstantools::use_rstan()` generates a `stanmodels` object/list; `get_stan_model()` becomes a simple lookup into that precompiled list.

### Phase 3 — Replace sampling, draws, diagnostics (parametric engine)
Refactor these functions to operate on `rstan` objects:

Affected functions:
- `R/bayesian-impute.R`
  - `bayesian_impute_single()`
  - `extract_posterior_distributions()`
  - `compute_diagnostics()`
  - `compute_waic_from_stan()`
- `R/utils.R`
  - `extract_mcmc_summary()`

Key API mappings (cmdstanr -> rstan):

1) Sampling
- Current: `stan_model$sample(iter_warmup, iter_sampling, chains, adapt_delta, max_treedepth, ...)`
- Target: `rstan::sampling(stan_model, data = stan_data, iter = warmup + sampling, warmup = warmup, chains = chains, control = list(adapt_delta = ..., max_treedepth = ...))`

2) Draws extraction for plotting (keep `posterior` compatibility)
- Current: `posterior_samples <- stan_fit$draws(param_names)` (a `posterior` draws object)
- Target: standardize to storing a `posterior` draws object regardless of backend, e.g.
  - Convert `stanfit` to a `posterior::draws_*` type (then existing plotting helpers keep working).

3) `t_imputed` extraction (posterior imputations)
- Current: `stan_fit$draws("t_imputed", format = "array")`
- Target (example approach): use `rstan::extract(stan_fit, pars = "t_imputed", permuted = FALSE)` then reshape to the same “total_draws x n_censored” matrix expected by downstream code.

4) Diagnostics
- Current: `diagnostic_summary()`, plus `summary(..., variables=...)` for Rhat/ESS
- Target:
  - Parameter summary: `rstan::summary(stan_fit, pars = param_names)$summary` (contains `Rhat`, `n_eff`)
  - Divergences / treedepth: `rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)` then count `divergent__` and `treedepth__` hits

5) WAIC
- Current: relies on Stan output variables `ll_obs` / `ll_cens` extracted via `stan_fit$draws(...)`
- Target:
  - Keep the same Stan model outputs (`ll_obs`, `ll_cens`) but replace the extraction with `rstan` extraction/conversion.
  - Alternatively (optional): adopt `loo` (CRAN) and compute WAIC/LOO from a single `log_lik` matrix if you adjust the Stan programs to emit `log_lik[n]`.

### Phase 4 — Tests + CI adjustments
Update tests to align with rstan:
- Replace `skip_if_not_installed("cmdstanr")` with `skip_if_not_installed("rstan")` *only if* rstan is not a hard dependency.
- If `rstan` becomes `Imports`, tests can assume it is present, but keep iteration counts small and `skip_on_cran()` for sampling-based tests.

Files to revisit:
- `tests/testthat/test-stan-models.R`
- `tests/testthat/test-export.R`
- `tests/testthat/test-complete.R`

CI considerations:
- GitHub Actions will need a setup that installs `rstan` (binaries are available on most runners but compilation can still be slow).
- Keep Stan sampling tests minimal, and prefer `skip_on_cran()` for anything that compiles/samples.

## File-by-file change map (what will change)

### Core changes
- `R/stan-models.R`: rewrite for `rstan` / `rstantools`
- `R/bayesian-impute.R`: replace sampling + draws extraction + diagnostics + WAIC extraction
- `R/utils.R`: remove `check_cmdstan()`, update `extract_mcmc_summary()`, update any “required packages” checks

### Metadata/docs
- `DESCRIPTION`: remove `cmdstanr`/`Additional_repositories`, add `rstan` (+ `LinkingTo` etc if A2)
- `README.md`: update installation/toolchain notes
- `man/check_cmdstan.Rd`: remove or replace (if the function is removed/renamed)

### Packaging/build (only for A2)
- `src/Makevars`, `src/Makevars.win`: merge rstantools requirements with existing Fortran build settings.
  - This is the riskiest integration point because the package already contains substantial compiled Fortran/C code.

## Open questions to resolve early

1) Do you require “install-time compiled Stan models” (A2), or is “compile-on-first-use” acceptable (A1)?
2) Should Stan-based parametric imputation remain the default engine in `impute()` (currently it is when model is not nonparametric)?
3) Which output type should `posterior_samples` be standardized to?
   - Keeping it as a `posterior::draws_*` object is ideal because plotting already expects it.

## Rough effort estimate (order-of-magnitude)

- **A1 (rstan runtime compilation):** ~1–2 weeks to implement + stabilize across platforms.
- **A2 (rstantools install-time compilation):** add ~1–3 weeks depending on how smoothly `src/` build integration goes (especially Windows/macOS).

The main schedule risk is platform build issues when combining rstantools-generated C++ compilation with the existing Fortran codebase in `src/`.

