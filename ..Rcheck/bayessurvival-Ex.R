pkgname <- "bayessurvival"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bayessurvival')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("apply_halabi_singh_censoring")
### * apply_halabi_singh_censoring

flush(stderr()); flush(stdout())

### Name: apply_halabi_singh_censoring
### Title: Apply Halabi-Singh Censoring Method
### Aliases: apply_halabi_singh_censoring

### ** Examples

## Not run: 
##D # Generate Weibull survival times
##D true_times <- simulate_survival_data(100, "weibull", list(shape = 2, scale = 4))
##D 
##D # Apply 30% exponential censoring 
##D censored_data <- apply_halabi_singh_censoring(true_times, censoring_rate = 0.3)
##D 
##D # Apply uniform censoring (for comparison)
##D censored_uniform <- apply_halabi_singh_censoring(
##D   true_times, censoring_rate = 0.3, censoring_distribution = "uniform"
##D )
##D 
##D # Check actual censoring rates
##D mean(censored_data$status == 0)
## End(Not run)



cleanEx()
nameEx("bayes_np_impute")
### * bayes_np_impute

flush(stderr()); flush(stdout())

### Name: bayes_np_impute
### Title: bayes_np_impute: Non-parametric Bayesian imputation for
###   right-censored survival data using a _Linear Dependent
###   Dirichlet-Process_ model (originally from *DPpackage*).
### Aliases: bayes_np_impute

### ** Examples

set.seed(1)
n  <- 50
t  <- rexp(n, 0.05)               # true times
c  <- rexp(n, 0.02)               # censor times
time   <- pmin(t, c)
status <- as.integer(t <= c)
df <- data.frame(time, status)

# Default: alpha ~ Gamma(10, 1) (learnable)
res <- bayes_np_impute(df, "time", "status", n_imputations = 10, verbose = FALSE,
                       mcmc = list(nburn = 100, nsave = 200, nskip = 2, ndisplay = 100))

# View results
print(res)
plot(res)

# Access completed datasets
head(res$imputed_datasets[[1]])

# Custom prior with fixed alpha = 1
custom_prior <- list(a0 = -1, b0 = -1, nu = 4, m0 = 0, 
                     S0 = matrix(100, 1, 1), psiinv = matrix(1, 1, 1),
                     tau1 = 6.01, taus1 = 6.01, taus2 = 2.01)
res_fixed <- bayes_np_impute(df, "time", "status", n_imputations = 5, 
                             prior = custom_prior, verbose = FALSE,
                             mcmc = list(nburn = 100, nsave = 200, nskip = 2, ndisplay = 100))




cleanEx()
nameEx("bayesian_impute")
### * bayesian_impute

flush(stderr()); flush(stdout())

### Name: bayesian_impute
### Title: Bayesian Imputation for Censored Survival Data
### Aliases: bayesian_impute

### ** Examples

## Not run: 
##D # Load example data
##D data(lung, package = "survival")
##D lung_clean <- na.omit(lung[, c("time", "status")])
##D 
##D # Basic imputation (single group)
##D result <- bayesian_impute(
##D   data = lung_clean,
##D   time_col = "time",  
##D   status_col = "status",  
##D   n_imputations = 20
##D )
##D 
##D # Group comparison
##D result_groups <- bayesian_impute(
##D   data = lung_clean,
##D   time_col = "time",
##D   status_col = "status",
##D   groups = "sex",
##D   n_imputations = 20
##D )
##D 
##D # View results
##D print(result)
##D plot(result, type = "survival")
##D 
##D # Generate additional datasets 
##D more_datasets <- generate_complete_datasets(
##D   result$original_data, 
##D   result$time_col, 
##D   result$status_col,
##D   result$posterior_imputations, 
##D   50  # Generate 50 more datasets
##D )
## End(Not run)



cleanEx()
nameEx("complete.bayesian_imputation")
### * complete.bayesian_imputation

flush(stderr()); flush(stdout())

### Name: complete.bayesian_imputation
### Title: Extract Complete Datasets from Bayesian Imputation
### Aliases: complete.bayesian_imputation

### ** Examples

## Not run: 
##D # Basic imputation
##D imp <- bayesian_impute(lung, "time", "status", n_imputations = 5)
##D 
##D # Extract first dataset
##D dataset1 <- complete(imp, dataset = 1)
##D 
##D # Extract all datasets
##D all_datasets <- complete(imp)  # or complete(imp, "all")
##D 
##D # Extract specific datasets
##D some_datasets <- complete(imp, dataset = c(1, 3, 5))
##D 
##D # Get long format (all datasets stacked)
##D long_data <- complete(imp, format = "long")
##D 
##D # Include original data with missing values
##D with_original <- complete(imp, include.original = TRUE)
## End(Not run)




cleanEx()
nameEx("dp_np_LDDPsurvival")
### * dp_np_LDDPsurvival

flush(stderr()); flush(stdout())

### Name: dp_np_LDDPsurvival
### Title: Non-parametric DP Survival (Wrapper)
### Aliases: dp_np_LDDPsurvival
### Keywords: internal

### ** Examples

# Not run: see bayes_np_impute() for a full example.
# dp_np_LDDPsurvival(ymat ~ 1, zpred = matrix(1,1,1), ...)




cleanEx()
nameEx("event_probability")
### * event_probability

flush(stderr()); flush(stdout())

### Name: event_probability
### Title: Event probability by clinical time horizons
### Aliases: event_probability

### ** Examples

# event_probability(result, times = c(180, 365))



cleanEx()
nameEx("event_probability_groups")
### * event_probability_groups

flush(stderr()); flush(stdout())

### Name: event_probability_groups
### Title: Event probability by clinical time horizons (groups)
### Aliases: event_probability_groups

### ** Examples

# event_probability_groups(result_groups, times = c(180, 365))



cleanEx()
nameEx("export")
### * export

flush(stderr()); flush(stdout())

### Name: export
### Title: Export Imputed Datasets
### Aliases: export

### ** Examples

## Not run: 
##D # Load and impute data
##D data(lung, package = "survival")
##D result <- bayesian_impute(lung, n_imputations = 5)
##D 
##D # Export all datasets as CSV files
##D export(result, "lung_imputed", format = "csv")
##D 
##D # Export as RDS for future R analysis
##D export(result, "lung_imputed", format = "rds")
##D 
##D # Export only first dataset
##D export(result, "lung_imputed", format = "csv", datasets = "first")
##D 
##D # For group analysis
##D result_groups <- bayesian_impute(data, groups = "treatment")
##D export(result_groups, "group_imputed", format = "csv", groups = "combined")
## End(Not run)




cleanEx()
nameEx("generate_complete_datasets")
### * generate_complete_datasets

flush(stderr()); flush(stdout())

### Name: generate_complete_datasets
### Title: Generate Complete Datasets from Posterior Distributions
### Aliases: generate_complete_datasets

### ** Examples

## Not run: 
##D # After running bayesian_impute
##D result <- bayesian_impute(data, "time", "status")
##D 
##D # Generate 50 additional datasets
##D more_datasets <- generate_complete_datasets(
##D   result$original_data,
##D   result$time_col,
##D   result$status_col,
##D   result$posterior_imputations,
##D   50
##D )
## End(Not run)



cleanEx()
nameEx("prepare_survival_data")
### * prepare_survival_data

flush(stderr()); flush(stdout())

### Name: prepare_survival_data
### Title: Prepare Survival Data for Bayesian Imputation
### Aliases: prepare_survival_data

### ** Examples

## Not run: 
##D # Auto-detection (recommended for beginners)
##D library(survival)
##D lung_data <- prepare_survival_data(lung)
##D 
##D # Manual specification  
##D my_data <- prepare_survival_data(mydata, time = "duration", status = "died")
##D 
##D # Then use directly with bayesian_impute
##D result <- bayesian_impute(lung_data, n_imputations = 20)
## End(Not run)




cleanEx()
nameEx("quick_export")
### * quick_export

flush(stderr()); flush(stdout())

### Name: quick_export
### Title: Quick Export Function
### Aliases: quick_export

### ** Examples

## Not run: 
##D result <- bayesian_impute(data, n_imputations = 5)
##D quick_export(result, "my_imputed_data")
## End(Not run)




cleanEx()
nameEx("simulate_survival_data")
### * simulate_survival_data

flush(stderr()); flush(stdout())

### Name: simulate_survival_data
### Title: Generate Survival Data from Known Distribution
### Aliases: simulate_survival_data

### ** Examples

## Not run: 
##D # Generate Weibull data 
##D times <- simulate_survival_data(100, "weibull", list(shape = 2, scale = 4))
##D 
##D # Generate exponential data
##D times <- simulate_survival_data(100, "exponential", list(rate = 0.5))
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
