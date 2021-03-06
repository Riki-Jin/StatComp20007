% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp20007R.R
\name{gaussnonr}
\alias{gaussnonr}
\title{Non-normal random numbers generator with desired Mardia skewness and kurtosis.}
\usage{
gaussnonr(p, n, ms, mk, sig, Sigma, initial = NULL)
}
\arguments{
\item{p}{Dimension}

\item{n}{Sample size}

\item{ms}{Desired Mardia skewness}

\item{mk}{Desired Mardia kurtosis}

\item{sig}{p*p multivariate normal covariance matrix}

\item{Sigma}{p*p desired result covariance matrix}

\item{initial}{A vector with 3 numbers for initial polynominal coefficients' (b,c,d). The default setting is (0.9,0.4,0).}
}
\value{
A data matrix (non-normal random numbers generated)
}
\description{
Non-normal random numbers generator with desired Mardia skewness and kurtosis.
}
\examples{
\dontrun{
set.seed(1)
p <- 2
n <- 10000
ms <- 3
mk <- 20
sig <- matrix(c(1,0.7,0.7,1),2,2)
Sigma <- matrix(c(1,0.5,0.5,1),2,2)
}
}
