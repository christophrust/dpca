#' S3 methods for s3 classes "dpca" and "spca"

#' Print method for object of class \code{dpca}.
#' @param x An object of type \code{dpca}.
#' @param ... Further pass-through arguments.
#' @export
print.dpca <- function(x, ...) {
  cat("\nDynamic principal component estimation object:\n")
  cat(sprintf("  Number of selected dynamic components: %s\n", x$HL_select$q))
}

#' Summary method for object of class \code{dpca}.
#' @param object An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @export
summary.dpca <- function(object, ...) print.dpca(object)

#' Print method for object of class \code{spca}.
#' @param x An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @export
print.spca <- function(x, ...) {
  cat("\nStatic principal component estimation object:\n")
  cat(sprintf("  Number of selected static components: %s\n", x$HL_select$r))
}

#' Summary method for object of class \code{spca}.
#' @param object An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @export
summary.spca <- function(object, ...) print.spca(object)

#' Plot q-selection graphcs for a \code{dpca} object.
#' @param x An object of type \code{spca}.
#' @param ... Further pass-through arguments.
#' @importFrom graphics axis mtext par
#' @export
plot.dpca <- function(x, ...) {
  if (is.null(x$HL_select)) {
    warning("No data-driven selection of number of dynamic factors for passed object. Nothing to plot!")
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
    type = "n", ylab = "Number of selected factors",
    xlab = "c", lwd = 2
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
