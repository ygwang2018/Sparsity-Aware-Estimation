\name{gdat}
\alias{gdat}
\alias{simu}
\alias{eval_methods}
\title{Simulation Utilities}
\description{Generate simulation data and run the simulation framework used in the manuscript.}
\usage{
gdat(n, K, p, rho, beta0, treat_prob = 0.3,
  balance = TRUE, hetero_error = TRUE)
}
\arguments{
  \item{n}{Total sample size.}
  \item{K}{Number of groups.}
  \item{p}{Number of covariates.}
  \item{rho}{AR(1)-type covariate correlation parameter.}
  \item{beta0}{Coefficient matrix in the orientation expected by \code{gdat}.}
  \item{treat_prob}{Marginal treatment probability.}
  \item{balance}{Whether group sizes are balanced.}
  \item{hetero_error}{Whether to use group-level heteroscedastic errors.}
}
\value{\code{gdat} returns a data frame with Y, A, G, and covariates, with true nuisance quantities stored as attributes.}
