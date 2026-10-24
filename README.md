# SparsityAwareCATE

R code for the manuscript **“Sparsity- and Homogeneity-Aware Estimation of Group-Specific Conditional Average Treatment Effects”** by Xin Ye and You-Gan Wang.

`SparsityAwareCATE` implements a penalized R-learning estimator that simultaneously:

- selects group-specific treatment-effect modifiers;
- identifies coefficient blocks that are homogeneous across all study groups;
- uses cross-fitted nuisance-function estimates;
- optionally uses group-specific loss weights to improve efficiency under heteroscedasticity;
- fits the structured penalty by ADMM; and
- provides oracle-style post-selection inference conditional on the selected sparsity and homogeneity structure.

## Installation

Install the package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("ygwang2018/Sparsity-Aware-Estimation")
```

Then load it with:

```r
library(SparsityAwareCATE)
```

The main package dependencies are `MASS`, `glmnet`, `grpreg`, `ipred`, and `randomForest`.

## Basic workflow

The input data frame should contain columns `Y`, `A`, and `G`, followed by baseline covariates. `A` is a binary treatment indicator and `G` is an integer group label from 1 to `K`.

```r
library(SparsityAwareCATE)

set.seed(2026)
K <- 4
p <- 20
n <- 2400

b1 <- c(-0.6, -0.6, 0.6, 0.6, 1.0, 0, 0, 0, rep(0, p - 8))
b2 <- c(-0.6, -0.6, 0.6, 0.6, 0, -0.8, 0.8, 0, rep(0, p - 8))
b3 <- c(-0.6, -0.6, 0.6, 0.6, 0, -0.8, 0.8, 0, rep(0, p - 8))
b4 <- c(-0.6, -0.6, 0.6, 0.6, 0, 0, 0, 1.0, rep(0, p - 8))
beta0 <- cbind(b1, b2, b3, b4)

dat <- gdat(
  n = n, K = K, p = p, rho = 0.5,
  beta0 = beta0,
  treat_prob = 0.5,
  balance = TRUE
)

# Cross-fitted nuisance estimation using bagging
mydat <- prepare_dat(
  dat,
  K = K,
  nfold = 10,
  outcome_method = "bagging",
  propensity_method = "bagging",
  outcome_args = list(nbagg = 50),
  propensity_args = list(nbagg = 50),
  seed = 2026,
  trim_e = 1e-4,
  use_group_for_nuisance = TRUE,
  add_group_interactions = FALSE,
  calibrate_propensity = TRUE,
  use_dc_screening = FALSE
)

# Initial group-specific estimator
init <- method_loc_lasso_ols_refit(mydat)

# Unweighted proposed estimator
fit <- method_prop_admm(
  mydat,
  beta_init = init$beta,
  wloss = rep(1, mydat$n)
)

fit$beta
fit$lambda1
fit$lambda2
```

## Optional loss weighting

Cross-fitted group-specific inverse-variance weights can be constructed by:

```r
wf <- make_wloss_cf_local(
  mydat,
  nfold = 10,
  seed = 2026
)

fit_w <- method_prop_admm(
  mydat,
  beta_init = init$beta,
  wloss = wf$wloss
)
```

These weights are intended as an efficiency-motivated practical option under group-level heteroscedasticity. The formal asymptotic theory in the manuscript treats the loss weights as bounded functions of the observed covariates and group.

## Oracle-style post-selection inference

The package provides structured sandwich inference conditional on a selected sparsity/homogeneity structure:

```r
structure_hat <- get_structure_from_beta(fit$beta)

# Contrast for coefficient 1 in group 1
contrast <- make_beta1_contrast(j = 1, K = K, p = p)

inf <- infer_contrast_structured(
  mydat = mydat,
  a = contrast,
  structure = structure_hat,
  beta_input = fit$beta,
  refit = TRUE,
  level = 0.95
)

inf[c("psi_hat", "se", "ci_low", "ci_high")]
```

The resulting intervals are **oracle-style post-selection intervals**: they condition on the selected structure and do not fully account for finite-sample uncertainty from variable or homogeneity selection.

## Comparison estimators

The repository also contains implementations of the comparison estimators used in the manuscript:

```r
method_loc_lasso(mydat)
method_pooled_lasso(mydat)
method_glasso(mydat)
method_hlasso(mydat)
method_R2lasso(mydat)
```

## Reproducing the manuscript simulations

The original code supplied for the simulation and inference studies is retained under:

```text
inst/reproducibility/original_code/estimation/
inst/reproducibility/original_code/inference/
inst/reproducibility/original_code/compare_nuisance/
```

After installation, locate these files with:

```r
system.file("reproducibility", "original_code", package = "SparsityAwareCATE")
```

The principal simulation driver is `simu()`. The `compare_nuisance/` directory contains the nuisance-function comparison scripts, including the random-forest specification with `ntree = 150` used in the revised manuscript. The supplied V2 scripts are retained unchanged under `inst/reproducibility/original_code/` so that the research code remains auditable.

## Data

The patient-level IST-3 data used in the manuscript are publicly available through Edinburgh DataShare:

https://doi.org/10.7488/ds/1350

The IST-3 patient-level data are **not redistributed in this repository**.

## Code availability

Repository: https://github.com/ygwang2018/Sparsity-Aware-Estimation

## Citation

Please cite the associated manuscript:

> Ye, X. and Wang, Y.-G. (2026). *Sparsity- and Homogeneity-Aware Estimation of Group-Specific Conditional Average Treatment Effects*. Manuscript.

## Authors

- Xin Ye
- You-Gan Wang (corresponding author), School of Statistics and Data Science, Guangdong University of Finance and Economics. Email: you-gan.wang@uq.edu.au

## License

GPL-3.
