% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recursive_filter.R
\name{recursive_filter}
\alias{recursive_filter}
\title{Multivariate recursive filter}
\usage{
recursive_filter(x, f)
}
\arguments{
\item{x}{A multiple time series supplied as a n by T matrix}

\item{f}{Filter coefficient matrices supplied as n by n by p array,
where p is the maximum lag order.}
}
\description{
Multivariate recursive filter
}
