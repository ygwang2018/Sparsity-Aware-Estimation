\name{infer_contrast_structured}
\alias{infer_contrast_structured}
\alias{get_structure_from_beta}
\alias{make_beta1_contrast}
\title{Oracle-Style Structured Post-Selection Inference}
\description{Construct sandwich inference for contrasts conditional on a selected sparsity and homogeneity structure.}
\usage{
infer_contrast_structured(mydat, a, structure, beta_input = NULL,
  beta0 = NULL, wloss = NULL, refit = TRUE, ridge = 1e-8,
  level = 0.95)

get_structure_from_beta(beta_hat, active_tol = 1e-6, homo_tol = 1e-6)

make_beta1_contrast(j, K, p)
}
\arguments{
  \item{mydat}{Prepared data returned by \code{prepare_dat}.}
  \item{a}{Contrast vector or K by p contrast matrix.}
  \item{structure}{List containing active and homogeneous structure indicators.}
  \item{beta_input}{Optional fitted coefficient matrix.}
  \item{beta0}{Optional true coefficient matrix for simulation diagnostics.}
  \item{wloss}{Optional loss weights.}
  \item{refit}{Whether to refit the selected oracle model before inference.}
  \item{ridge}{Small numerical ridge used if inversion is unstable.}
  \item{level}{Confidence level.}
  \item{beta_hat}{Estimated K by p coefficient matrix.}
  \item{active_tol, homo_tol}{Numerical classification tolerances.}
  \item{j}{Covariate index.}
  \item{K}{Number of groups.}
  \item{p}{Number of covariates.}
}
\value{For inference, a list containing the contrast estimate, standard error, variance, and confidence interval.}
\note{Inference is oracle-style and conditions on the selected sparsity and homogeneity structure.}
