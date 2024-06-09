
#' Stacked bar plot
#' 
#' Produces stacked bar plots using base graphics.
#' 
#' @param x matrix of deconvolution results with samples in rows and cell
#'   subclasses or groups in columns.
#' @param percent Logical whether to scale the matrix rows as percentage.
#' @param order_col Numeric value for which column to sort. If a vector is
#'   supplied, these columns are averaged first using `rowMeans()`.
#' @param scheme Vector of colours. If not supplied, the default scheme uses
#'   `scales::hue_pal()`.
#' @param ... Optional arguments passed to [graphics::barplot()].
#' @returns No return value. Plots a stacked barchart using base graphics.
#' @importFrom scales hue_pal
#' @importFrom graphics barplot par
#' @export

stack_plot <- function(x, percent = FALSE, order_col = 1, scheme = NULL,
                       ...) {
  if (is.null(scheme)) {
    scheme <- hue_pal(h = c(0, 270))(ncol(x))
  }
  if (percent) {
    rs <- rowSums(x)
    x <- x / rs * 100
    ord <- if (length(order_col) == 1) {
      order(x[, order_col])
    } else {
      order(rowMeans(x[, order_col, drop = FALSE]))
    }
  } else {
    ord <- order(rowSums(x))
  }
  op <- par(mar = c(8, 4, 1.5, 1.5))
  on.exit(par(op))
  barplot(t(x[ord,]), las = 2, col = scheme,
          cex.names = 0.7,
          tcl = -0.3, mgp = c(2.2, 0.5, 0))
}


