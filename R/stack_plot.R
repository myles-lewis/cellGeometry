
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
#' @param cex.names Character expansion controlling bar names font size.
#' @param ... Optional arguments passed to [graphics::barplot()].
#' @returns No return value. Plots a stacked barchart using base graphics.
#' @importFrom scales hue_pal
#' @importFrom graphics axis barplot par strwidth
#' @export

stack_plot <- function(x, percent = FALSE, order_col = 1, scheme = NULL,
                       cex.names = 0.7, ...) {
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
  strw <- max(strwidth(rownames(x), units = "inches", cex = cex.names), na.rm = TRUE)
  mar1 <- strw / par("csi") +1.5
  op <- par(mar = c(mar1, 4, 1.5, 1.5))
  on.exit(par(op))
  barplot(t(x[ord,]), las = 2, col = scheme,
          cex.names = cex.names,
          tcl = -0.3, mgp = c(2.2, 0.5, 0), ...)
  # extend axis line
  yrange <- par("usr")[3:4]
  pos <- par("usr")[1]
  axis(2, at = c(0, yrange[2]), pos = pos, labels = FALSE, lwd.ticks = 0)
}


