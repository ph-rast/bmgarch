% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{.get_stan_summary}
\alias{.get_stan_summary}
\title{Get stan summaries.}
\usage{
.get_stan_summary(model_fit, params, CrI, weights = NULL, sampling_algorithm)
}
\arguments{
\item{model_fit}{stanfit object or list of stanfit objects.}

\item{params}{Character vector. Names of params to pull from stan summary.}

\item{CrI}{Numeric vector (length 2).}

\item{weights}{Numeric vector. Weights for each model in model_fit, if list.}

\item{sampling_algorithm}{Character vector for sampling method.}
}
\value{
Stan summary for parameters. Columns: mean, sd, mdn, and CrIs.
}
\description{
Get stan summaries.
}
\author{
Stephen R. Martin, Philippe Rast
}
\keyword{internal}
