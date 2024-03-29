% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SparseICA-package.R
\name{est.M.ols}
\alias{est.M.ols}
\title{Estimate mixing matrix from estimates of components}
\usage{
est.M.ols(sData, xData, intercept = TRUE)
}
\arguments{
\item{sData}{S Dimension: P x Q}

\item{xData}{X Dimension: P x T}

\item{intercept}{default = TRUE}
}
\value{
a mixing matrix M, dimension Q x T.
}
\description{
Estimate mixing matrix from estimates of components
}
