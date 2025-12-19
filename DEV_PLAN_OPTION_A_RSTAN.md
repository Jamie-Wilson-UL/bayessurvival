# Development Plan: Option A (CRAN) — migrate from `cmdstanr` to `rstan` / `rstantools`

This plan breaks the Option A migration into small, verifiable tasks in a safe order. It assumes **`rstan` becomes a hard dependency** (`Imports`) because Stan-based parametric imputation is core functionality.

## Progress tracker

Last updated: 2025-12-17

### Phase 0 — Baseline and CRAN hygiene

- [x] 0.1 Capture baseline check output (before changes)
  - Built tarball: `bayessurvival_0.1.0.9000.tar.gz`
  - `R CMD check --as-cran` result: failed due to (a) no internet for incoming/URL checks and (b) missing suggested packages (`writexl`, `covr`, `DT`, `flexsurv`), see `bayessurvival.Rcheck/00check.log`
- [x] 0.2 Remove bundled binaries/build artifacts (CRAN blocker)
  - Removed committed compiled artifacts from `src/` and removed the Stan executable `inst/stan/weibull_imputation`
- [x] 0.3 Remove accidental/stray files that can enter the tarball
  - Removed `tests/.DS_Store` and `lung_imputed_dataset_1.csv`
  - Hardened ignores in `.Rbuildignore` and `.gitignore`
  - Re-check note: `R CMD check --no-manual` passes (with `_R_CHECK_FORCE_SUGGESTS_=false`); Makevars portability warnings have been resolved

### Phase 1 — Metadata switch: `cmdstanr` → `rstan`

- [x] 1.1 Update `DESCRIPTION` (CRAN-facing)
  - Removed `Additional_repositories`
  - Removed `cmdstanr` and CmdStan-specific `SystemRequirements`
  - Added `rstan` to `Imports`
- [x] 1.2 Update README/docs to reflect `rstan`
  - Updated `README.md` to describe the `rstan` backend and point to `DEV_PLAN_OPTION_A_RSTAN.md`

### Phase 2 — Backend refactor

- [x] 2.1 Implement internal `rstan` model compilation/caching
  - Updated `R/stan-models.R` to compile via `rstan::stan_model()` and cache in `.bayessurvival_env`
- [x] 2.2 Remove/replace CmdStan checks
  - Replaced `check_cmdstan()` with `check_rstan()` and removed `man/check_cmdstan.Rd`

### Phase 3 — Implement A1 (rstan compile-on-first-use)

- [x] 3.1 Replace sampling call in `bayesian_impute_single()`
  - Updated `R/bayesian-impute.R` to use `rstan::sampling()`
- [x] 3.2 Replace draws extraction / posterior imputations / WAIC extraction
  - Updated `R/bayesian-impute.R` to use `posterior::as_draws_df()` and `rstan::extract()` for `t_imputed`, `ll_obs`, `ll_cens`
- [x] 3.3 Replace diagnostics extraction (`compute_diagnostics()`)
  - Updated `R/bayesian-impute.R` to use `posterior` for Rhat/ESS and `rstan::get_sampler_params()` for divergences/treedepth
- [x] 3.4 Update MCMC summary (`extract_mcmc_summary()`)
  - Updated `R/utils.R` to work with `stanfit` objects via `posterior`

### Phase 4 — Tests/CI alignment

- [x] 4.1 Update tests away from `cmdstanr`/CmdStan
  - Updated `tests/` to remove CmdStan assumptions and stop skipping mock tests
- [x] 4.2 Verify checks (local)
  - `R CMD build` + `R CMD check --no-manual` passes (with `_R_CHECK_FORCE_SUGGESTS_=false`)
  - Remaining NOTEs: suggested packages unavailable, empty `NEWS.md`, Fortran I/O note from vendored DPpackage code

## Success criteria

- `DESCRIPTION` contains no `Additional_repositories` and no `cmdstanr`.
- Parametric `bayesian_impute()` runs using `rstan` with the existing `.stan` programs in `inst/stan/`.
- `R CMD build` + `R CMD check --as-cran` on the **built tarball** has **0 ERROR, 0 WARNING** (only explainable NOTEs).
- Tests pass locally; CI is green on at least macOS + Windows + Linux (release).

## Phase 0 — Baseline and CRAN hygiene (blockers first)

### 0.1 Capture baseline check output (before changes)
- Run: `R CMD build .`
- Run: `R CMD check --as-cran <tarball>`
- Save notes: current WARN/NOTE list to track regressions.

Gate:
- You have a baseline `00check.log` to compare against later.

### 0.2 Remove bundled binaries/build artifacts (CRAN blocker)
- Delete any committed compilation outputs:
  - `src/*.o`, `src/*.so` (and any platform variants)
  - Stan executables under `inst/stan/` (e.g., `inst/stan/weibull_imputation`)
- Ensure these are not committed again (gitignore + build practices).

Gate:
- `R CMD build` tarball does **not** include object/shared-library artifacts in `src/`.

### 0.3 Remove accidental/stray files that can enter the tarball
- Remove junk files like `tests/.DS_Store`.
- Remove or relocate any root-level artifacts not meant for CRAN (e.g., exported CSV outputs).

Gate:
- `tar -tf <tarball> | rg -n "(DS_Store|\\.csv$|\\.Rcheck)"` only shows intentional package files.

## Phase 1 — Metadata switch: `cmdstanr` → `rstan`

### 1.1 Update `DESCRIPTION` (CRAN-facing)
- Remove:
  - `Additional_repositories`
  - `cmdstanr` from `Suggests`
  - CmdStan-specific `SystemRequirements` wording
- Add:
  - `Imports: rstan` (and any other required CRAN packages you adopt during the refactor)
- Optional (later, if you adopt `rstantools`): add required `LinkingTo`/`Imports` entries per scaffold.

Gate:
- `R CMD check --as-cran` does not report non-CRAN repositories or missing dependencies.

### 1.2 Update README/docs to reflect `rstan`
- Replace CmdStan installation instructions with `rstan` toolchain guidance.
- Ensure examples do not force long compilation or sampling during checks.

Gate:
- No documentation instructs installing from non-CRAN repos for core use.

## Phase 2 — Backend refactor (create a stable internal interface)

### 2.1 Define internal “Stan backend” helpers (rstan implementation)
Create/standardize internal helpers so the rest of the code doesn’t care about the backend:
- `get_stan_model(distribution)` → returns compiled model object
- `run_stan_sampling(model, data, mcmc_options)` → returns fit object
- `extract_draws_matrix(fit, var)` and/or `extract_draws_array(fit, var)`
- `extract_param_summary(fit, vars)` (Rhat/ESS)
- `extract_sampler_diagnostics(fit)` (divergences, treedepth hits)

Gate:
- A tiny “smoke call” can compile and sample one model on your machine.

### 2.2 Remove/replace CmdStan checks
- Remove `check_cmdstan()` and its Rd entry (or replace with `check_rstan()`).
- Update any references in code/tests/docs.

Gate:
- `R CMD check` has no orphaned man pages (no dangling `man/check_cmdstan.Rd` issues).

## Phase 3 — Implement A1 (rstan compile-on-first-use) end-to-end

### 3.1 Replace model compilation/caching (`R/stan-models.R`)
- Replace `cmdstanr::cmdstan_model()` with `rstan::stan_model(file=...)` using:
  - `system.file("stan", "<model>.stan", package = "bayessurvival")`
- Keep caching in `.bayessurvival_env` (same approach as now).

Gate:
- `get_stan_model("weibull")` returns a compiled `stanmodel` without errors.

### 3.2 Replace sampling call in `bayesian_impute_single()`
Map options:
- Current: `iter_warmup`, `iter_sampling`, `chains`, `adapt_delta`, `max_treedepth`
- Target:
  - `rstan::sampling(model, data = stan_data, chains = ..., iter = warmup + sampling, warmup = warmup, control = list(adapt_delta=..., max_treedepth=...))`

Gate:
- A tiny run (e.g., warmup=50, sampling=50, chains=2) finishes and returns a fit.

### 3.3 Replace draws extraction
Update these call sites:
- `posterior_samples <- stan_fit$draws(param_names)`
- `extract_posterior_distributions()` uses `stan_fit$draws("t_imputed", format="array")`
- `compute_waic_from_stan()` uses `stan_fit$draws(...)`

Target behavior to preserve:
- `posterior_imputations`: a matrix `[total_draws x n_censored]`
- `posterior_samples`: store as a `posterior::draws_*` object if feasible so plotting continues to work as-is.

Gate:
- `plot_trace_plots()` and `plot_pairs()` still work on a result object.

### 3.4 Replace diagnostics extraction (`compute_diagnostics()`)
- Replace cmdstanr diagnostics with `rstan` equivalents:
  - `rstan::summary(fit, pars=...)$summary` (Rhat / n_eff)
  - `rstan::get_sampler_params(fit, inc_warmup = FALSE)` (divergences, treedepth hits)

Gate:
- `result$diagnostics$convergence_ok` evaluates without error and fields remain populated.

### 3.5 Replace WAIC extraction (`compute_waic_from_stan()`)
- Keep your current WAIC computation logic.
- Change only how `ll_obs`/`ll_cens` are extracted from the fit.

Gate:
- If Stan outputs `ll_obs`/`ll_cens`, WAIC computes; if not present, behavior is unchanged (error/NULL per your design).

## Phase 4 — Tests and CI alignment

### 4.1 Update test skips and assumptions
- Replace `skip_if_not_installed("cmdstanr")` with `skip_if_not_installed("rstan")` only where needed.
- Keep `skip_on_cran()` for sampling/compile tests to avoid timeouts.
- Since `rstan` is an `Imports`, prefer tests that assume availability (but still avoid heavy compilation in CRAN checks).

Gate:
- `devtools::test()` passes locally from a clean install.

### 4.2 Add one minimal parametric integration test (recommended)
- Use a tiny dataset + small iterations to validate the full pipeline returns:
  - `posterior_imputations` shape
  - completed dataset generation
  - key fields present
- Mark with `skip_on_cran()` and keep runtime small.

Gate:
- Integration test is stable and doesn’t materially slow local runs.

### 4.3 Update GitHub Actions workflows
- Ensure Linux runners can install `rstan` system dependencies.
- Ensure Windows/macOS jobs remain green.

Gate:
- CI matrix passes `R CMD check` jobs.

## Phase 5 (optional) — Move from A1 to A2 (`rstantools`)

Only start once Phase 3–4 are stable.

### 5.1 Adopt rstantools scaffolding
- Integrate `rstantools` so models build at install time.
- Carefully merge build changes with existing Fortran/C compilation in `src/` (highest risk area).

Gate:
- Package installs from source on Windows/macOS/Linux without manual intervention.

## Suggested “stop points” for review

- After Phase 1: metadata/docs CRAN-friendly
- After Phase 3.2: rstan sampling works end-to-end
- After Phase 3.4: diagnostics match previous behavior
- After Phase 4: tests/CI stable; ready for pre-CRAN external checks
