\name{method_prop_admm}
\alias{method_prop_admm}
\alias{prepare_dat}
\alias{make_wloss_cf_local}
\alias{method_loc_lasso}
\alias{method_loc_lasso_ols_refit}
\alias{method_pooled_lasso}
\alias{method_glasso}
\alias{method_hlasso}
\alias{method_R2lasso}
\alias{method_oracle}
\title{Core Estimation Functions}
\description{Core estimation and preprocessing functions for group-specific CATE analysis.}
\usage{
prepare_dat(dat, K, nfold = 5, outcome_method = "lm",
  propensity_method = "glm", outcome_args = NULL,
  propensity_args = NULL, seed = 1, trim_e = 1e-4,
  use_group_for_nuisance = TRUE, add_group_interactions = TRUE,
  calibrate_propensity = TRUE, use_dc_screening = TRUE,
  dc_screen_n = 500, dc_tail_n = 30, dc_z = 3,
  dc_min_selected = 1, dc_seed = seed, dc_verbose = FALSE)

method_prop_admm(mydat, beta_init = NULL, wloss = NULL, gamma = 1,
  nlambda1 = 10, nlambda2 = 10, lambda.min.ratio = 1e-4,
  lambda1.max.multiplier = 10, lambda2.max.multiplier = 10,
  lambda1_seq = NULL, lambda2_seq = NULL, eps_weight = 1e-4,
  rho1 = 1, rho2 = 1, maxit = 1000, eps_pri = 1e-4,
  eps_dual = 1e-4, df_tol = 1e-6, verbose = FALSE,
  zero_init_tol = 0)
}
\arguments{
  \item{dat}{Data frame containing Y, A, G, and baseline covariates.}
  \item{K}{Number of groups.}
  \item{mydat}{Prepared data returned by \code{prepare_dat}.}
  \item{beta_init}{Initial K by p coefficient matrix.}
  \item{wloss}{Optional observation-level loss weights.}
  \item{nfold, outcome_method, propensity_method, outcome_args, propensity_args, seed, trim_e, use_group_for_nuisance, add_group_interactions, calibrate_propensity, use_dc_screening, dc_screen_n, dc_tail_n, dc_z, dc_min_selected, dc_seed, dc_verbose}{Arguments controlling nuisance estimation and optional screening.}
  \item{gamma, nlambda1, nlambda2, lambda.min.ratio, lambda1.max.multiplier, lambda2.max.multiplier, lambda1_seq, lambda2_seq, eps_weight, rho1, rho2, maxit, eps_pri, eps_dual, df_tol, verbose, zero_init_tol}{Arguments controlling adaptive weights, tuning grids, and the ADMM solver.}
}
\value{A list containing fitted coefficient estimates and method-specific diagnostics.}
\seealso{\code{infer_contrast_structured}}
