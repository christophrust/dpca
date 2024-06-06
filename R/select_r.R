#' Select number of static factors.
#'
#' Select the number of "static" factors for a static factor sequence
#' via the criteria of Bai and Ng (2002) but using the stability interval
#' method suggested by Hallin and Liska (2007).
#'
#' @param x Input data supplied either as a matrix (rows correspond to cross-sectional
#' units and columnts to observations in the time domain) or a multivariate object of
#' class \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}.
#'
#' @param crit Either of \code{"IC1"}, \code{"IC2"}, \code{"IC3"},
#' specifying which penalty to use. See Bai and Ng (2002) for
#' details on the criteria. Defaults to \code{"IC1"}
#'
#' @param penalty_scales A vector of penalty scales over which the stability
#' is evaluated. See Hallin and Liska (2007) for details.
#'
#' @param n_path Integer vector specifying which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure.
#'
#' @param max_r Integer, maximum number of components considered.
#'
#' @param ... Further arguments passed to internal
#' methods, currently without functionality.
#'
#' @return A list with the entries
#' \itemize{
#'   \item \code{evals}: the first \code{r} values are the eigenvalues
#'     of the covariance matrix
#'   \item \code{evecs}: the first \code{r} eigenvectors of the
#'     covariance matrix
#'   \item \code{unpenalized_ic_vals}: Unpenalized values of the selected
#'     information criterion.
#'   \item \code{sample_var_criterion}: \code{sample variance} of the
#'     selected \code{r} overall entries in \code{n_path} and
#'     \code{penalty_scales}.
#'   \item \code{info}: single integer indicating success/failure finding
#'     a stability interval. It can take the following values:
#'      0: everything went fine.
#'      1: no zero stability invervals were found.
#'      2: no stability was found, such that the penalty scale which globally
#'         minimizes the sample variance is chosen.
#'   \item \code{r}: the number of selected components.
#' }
#'
#' @references Hallin, M. and Liska, R. (2007). Determining the Number of
#' Factors in the General Dynamic Factor Model. Journal of the American
#' Statistical Association, 102 (478).
#'
#' @importFrom stats cov
#' @export
select_r <- function(
  x,
  crit = c("IC1", "IC2", "IC3"),
  penalty_scales = seq(0, 2, by = 0.01),
  n_path,
  max_r,
  ...
) {

  x <- if (is.ts(x) || "zoo" %in% class(x)) {
         t(x)
       } else if (is.matrix(x)) {
         x
       } else {
         stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
       }

  n <- nrow(x)
  t_len <- ncol(x)
  cx <- cov(t(x))

  if (missing(n_path)) {
    n_path <- floor(seq(n / 2, n, n / 20))
  }

  if (missing(max_r)) {
    max_r <- floor(sqrt(n))
  }

  ## penalties suggested by Bai & Ng (2002)
  penalties <- if (crit == "IC1") {
    (n_path + t_len) / (n_path * t_len) *
      log((n_path * t_len) / (n_path + t_len))
  } else if (crit == "IC2") {
    (n_path + t_len) / (n_path * t_len) * log(pmin(n_path, t_len))
  } else if (crit == "IC3") {
    log(pmin(n_path, t_len)) / (pmin(n_path, t_len))
  } else {
    stop("crit must be either of \"IC1\", \"IC2\", or \"IC3\"")
  }

  res <- .Call(
    "R_hl_select_q",
    cx,
    n_path,
    max_r,
    n,
    1L,
    2L,
    .Machine$double.eps,
    penalties,
    penalty_scales,
    PACKAGE = "dpca"
  )

  names(res)[names(res) == "q"] <- "r"
  res
}
