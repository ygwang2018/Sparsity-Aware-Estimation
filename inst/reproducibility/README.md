# Reproducibility code

This directory contains the authors' supplied simulation and inference scripts used during development of the manuscript.

- `original_code/estimation/`: estimation and simulation code.
- `original_code/inference/`: post-selection inference simulation code.

The installable package functions in `R/` are assembled from these supplied sources. Hyperparameters for nuisance learners can be passed explicitly to `prepare_dat()` through `outcome_args` and `propensity_args`.

For the nuisance-function comparison described in the manuscript, the reported settings are:

- bagging: `ipred`, `nbagg = 50`;
- random forest: `randomForest`, `ntree = 150`;
- lasso nuisance model: `glmnet`, `alpha = 1`, with five-fold internal cross-validation;
- nuisance estimation: ten-fold cross-fitting;
- propensity estimates truncated to `[1e-4, 1 - 1e-4]`.

Example:

```r
mydat_rf <- prepare_dat(
  dat,
  K = K,
  nfold = 10,
  outcome_method = "rf",
  propensity_method = "rf",
  outcome_args = list(ntree = 150),
  propensity_args = list(ntree = 150),
  trim_e = 1e-4,
  use_group_for_nuisance = TRUE,
  add_group_interactions = FALSE,
  calibrate_propensity = TRUE,
  use_dc_screening = FALSE
)
```
