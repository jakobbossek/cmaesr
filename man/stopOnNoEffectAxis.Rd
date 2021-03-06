% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stoppingConditions.R
\name{stopOnNoEffectAxis}
\alias{stopOnNoEffectAxis}
\title{Stopping condition: principal axis.}
\usage{
stopOnNoEffectAxis()
}
\value{
[\code{cma_stopping_condition}]
}
\description{
Stop if addition of 0.1 * sigma in a principal axis
direction does not change mean value.
}
\seealso{
Other stopping.conditions: \code{\link{stopOnCondCov}},
  \code{\link{stopOnMaxIters}},
  \code{\link{stopOnNoEffectCoord}},
  \code{\link{stopOnOptParam}},
  \code{\link{stopOnOptValue}},
  \code{\link{stopOnTimeBudget}}
}

