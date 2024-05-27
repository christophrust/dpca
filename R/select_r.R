#' Select number of static factors.
#'
#' Select the number of "static" factors for a static factor sequence
#' via the criteria of Bai and Ng (2002) but using the stability interval
#' method suggested by Hallin and Liska (2007).
#'
#' @param x A data frame of variables (each row one time observation)
#' or a data matrix (n by T).
#' @param crit Either of "IC1", "IC2", "IC3", specifying which penalty
#' to use. See Bai and Ng (2002) for details on the criteria.
#' Defaults to "IC1"
#' @param penalty_scales A vector of penalty scales over which the stability
#' is evaluated. See Hallin and Liska (2007) for details.
#' @param n_path Integer vector specifying which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure.
#' @param max_r Integer, maximum number of components considered.
#' @param ... Further arguments passed to internal
#' methods, currently without functionality.
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

  if (is.data.frame(x)) {

    n <- ncol(x)
    t_len <- nrow(x)
    cx <- cov(x)

  } else if (is.matrix(x)) {

    n <- nrow(x)
    t_len <- ncol(x)
    cx <- cov(t(x))

  } else {
    stop("x must be either a data frame, a data matrix (n by T)!")
  }

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
