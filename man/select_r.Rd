% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_r.R
\name{select_r}
\alias{select_r}
\title{Select number of static factors.}
\usage{
select_r(
  x,
  crit = c("IC1", "IC2", "IC3"),
  penalty_scales = seq(0, 2, by = 0.01),
  n_path,
  max_r,
  ...
)
}
\arguments{
\item{x}{Input data supplied either as a matrix (rows correspond to cross-sectional
units and columnts to observations in the time domain) or a multivariate object of
class \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}.}

\item{crit}{Either of \code{"IC1"}, \code{"IC2"}, \code{"IC3"},
specifying which penalty to use. See Bai and Ng (2002) for
details on the criteria. Defaults to \code{"IC1"}}

\item{penalty_scales}{A vector of penalty scales over which the stability
is evaluated. See Hallin and Liska (2007) for details.}

\item{n_path}{Integer vector specifying which (nested) subsets of the
cross section are used in the Hallin & Liska procedure.}

\item{max_r}{Integer, maximum number of components considered.}

\item{...}{Further arguments passed to internal
methods, currently without functionality.}
}
\value{
A list with the entries
\itemize{
  \item \code{evals}: the first \code{r} values are the eigenvalues
    of the covariance matrix
  \item \code{evecs}: the first \code{r} eigenvectors of the
    covariance matrix
  \item \code{unpenalized_ic_vals}: Unpenalized values of the selected
    information criterion.
  \item \code{sample_var_criterion}: \code{sample variance} of the
    selected \code{r} overall entries in \code{n_path} and
    \code{penalty_scales}.
  \item \code{info}: single integer indicating success/failure finding
    a stability interval. It can take the following values:
     0: Everything went fine.
     1: No zero stability invervals were found and the stability interval
        with lowest sample variance was chosen.
     2: No stability was found at all , such that the penalty scale which
        globally minimizes the sample variance is chosen.
  \item \code{r}: the number of selected components.
}
}
\description{
Select the number of "static" factors for a static factor sequence
via the criteria of Bai and Ng (2002) but using the stability interval
method suggested by Hallin and Liska (2007).
}
\references{
Hallin, M. and Liska, R. (2007). Determining the Number of
Factors in the General Dynamic Factor Model. Journal of the American
Statistical Association, 102 (478).
}
