#' S3 methods for s3 classes "dpca" and "spca"

#' Print method for object of class \code{dpca}.
#' @param x An object of type \code{dpca}.
#' @param ... Further pass-through arguments.
#' @export
print.dpca <- function(x, ...) {
  cat("\nDynamic principal component estimation\n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  cat(paste("Number of dynamic components:", dim(x$ndpc)[1]), "\n")
}

#' Summary method for object of class \code{dpca}.
#' @param object An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @export
summary.dpca <- function(object, ...) {
  ps <- object$eig$vectors
  lambda <- object$eig$values
  ndim <- dim(ps)[2]
  spec_chi <- vapply(seq_along(object$freqs), function(i) {
    crossprod(
      Conj(ps[, , i]) * lambda[, i], ps[, , i]
    )
  }, matrix(0i, nrow = ndim, ncol = ndim))

  z <- list()
  z$freqs <- object$freqs
  z$gamma <- apply(object$spectrum, c(1, 2), sum) / length(object$freqs)
  z$gamma_chi <- apply(spec_chi, c(1, 2), sum) / length(object$freqs)
  z$q <- dim(ps)[1]
  z$object <- object

  class(z) <- "summary.dpca"
  z
}

#' Summary method for object of class \code{summary.dpca}.
#'
#' @param x An object of type \code{summary.dpca}.
#' @param ... Further pass-through arguments.
#' @export
print.summary.dpca <- function(x, ...) {
  trace_chi <- Re(sum(diag(x$gamma_chi)))
  trace_x <- Re(sum(diag(x$gamma)))
  cat("\nDynamic principal component estimation summary\n\n")
  cat("\nCall:\n", paste(deparse(x$object$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  cat(paste("Number of dynamic components:", x$q, "\n"))
  cat(paste("Proportion of explained variance: ", round(trace_chi / trace_x, 2), "\n"))
}

#' Print method for object of class \code{spca}.
#' @param x An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @export
print.spca <- function(x, ...) {
  cat("\nStatic principal component estimation\n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  cat(paste("Number of dynamic components:", dim(x$factors)[1]), "\n")
}

#' Summary method for object of class \code{spca}.
#' @param object An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @export
summary.spca <- function(object, ...) print.spca(object)

#' Diagnostics plot for a \code{dpca} object.
#'
#' @param x An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @importFrom graphics axis mtext par
#' @export
plot.dpca <- function(x, ...) {
  if (is.null(x$HL_select)) {
    warning(paste(
      "No data-driven selection of number of dynamic factors for passed object.",
      "Nothing to plot!"
    ))
    return(invisible())
  }
  par(mar = c(5, 4, 4, 6))
  plot(
    x = x$HL_select$penalty_scales, y = x$HL_select$q_path,
    type = "n", col = "blue", ylab = "Number of selected factors",
    xlab = "c", lwd = 2, xlim = c(0, 1)
  )
  par(new = TRUE)
  plot(
    x = x$HL_select$penalty_scales, y = x$HL_select$sample_var, type = "l",
    axes = FALSE, bty = "n", xlab = "", ylab = "",
    col = "red", lwd = 1.5, xlim = c(0, 1)
  )
  axis(4)
  mtext(expression(paste("Sample variance of ", S[c])), side = 4, padj = 3)
  par(new = TRUE)
  plot(
    x = x$HL_select$penalty_scales, y = x$HL_select$q_path,
    axes = FALSE, bty = "n", xlab = "", ylab = "",
    type = "l", col = "blue", lwd = 3, xlim = c(0, 1)
  )
}

#' @export
plot.spca <- function(x, ...) {
  r_selection <- x$HL_select
  plot(
    x = r_selection$penalty_scales, y = r_selection$q_path,
    type = "n", ylab = "Number of factors",
    xlab = "c", lwd = 2,
    main = "Hallin & Liska (2007) factor selection summary"
  )
  par(new = TRUE)
  plot(
    x = r_selection$penalty_scales, y = r_selection$sample_var, type = "l",
    axes = FALSE, bty = "n", xlab = "", ylab = "",
    col = "red", lwd = 1.5
  )
  axis(4)
  mtext(expression(paste("Sample variance of ", S[c])), side = 4, padj = 3)
  par(new = TRUE)
  plot(
    x = r_selection$penalty_scales, y = r_selection$q_path,
    axes = FALSE, bty = "n", xlab = "", ylab = "",
    type = "l", col = "blue", lwd = 3
  )
}
