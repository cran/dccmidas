% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{riskmetrics_est}
\alias{riskmetrics_est}
\title{RiskMetrics model}
\usage{
riskmetrics_est(rt, lambda = 0.94)
}
\arguments{
\item{rt}{List of daily returns}

\item{lambda}{\strong{optional} Decay parameter. Default to 0.94}
}
\value{
A list with the \eqn{H_t} matrix, for each \eqn{t}.
}
\description{
Obtains the matrix H_t, under the RiskMetrics model.
}
\keyword{internal}
