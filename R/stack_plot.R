
#' Stacked bar plot
#' 
#' Produces stacked bar plots using base graphics or ggplot2 showing amounts of
#' cell subclasses in deconvoluted bulk samples.
#' 
#' @param x matrix of deconvolution results with samples in rows and cell
#'   subclasses or groups in columns. If a 'deconv' class object is supplied the
#'   deconvolution values for the cell subclasses are extracted and plotted.
#' @param percent Logical whether to scale the matrix rows as percentage.
#' @param order_col Numeric value for which column (cell subclass) to use to
#'   sort the bars - this only applies if `percent = TRUE`. If a vector of
#'   column indices is supplied, these columns are averaged first using
#'   `rowMeans()`. If `percent = FALSE`, then the default is to sort bars from
#'   low to high based on the row sums (i.e. total subclass cell amounts in each
#'   sample). Setting `order_col = 0` disables sorting of bars; in this case
#'   bars are shown in the original order of the rows of `x`.
#' @param scheme Vector of colours. If not supplied, the default scheme uses
#'   `scales::hue_pal()`.
#' @param cex.names Character expansion controlling bar names font size.
#' @param order_cells Character value specifying with cell types are ordered by
#'   abundance.
#' @param seriate Character value which enables ordering of samples using the
#'   `seriation` package. Any matrix based seriation methods can be used to
#'   order the samples. Recommended options include "CA", "BEA" or "BEA_TSP".
#' @param show_xticks Logical whether to show rownames as x axis labels.
#' @param legend_ncol Number of columns for ggplot2 legend. If set to `NULL`
#'   ggplot2 sets the column number automatically.
#' @param legend_position Position of ggplot2 legend
#' @param ... Optional arguments passed to [graphics::barplot()].
#' @returns The base graphics function has no return value. It plots a stacked
#'   barchart using base graphics. The ggplot2 version returns a ggplot2 object.
#' @importFrom scales hue_pal
#' @importFrom graphics axis barplot par strwidth
#' @export

stack_plot <- function(x, percent = FALSE, order_col = 1, scheme = NULL,
                       order_cells = c("none", "increase", "decrease"),
                       seriate = NULL,
                       cex.names = 0.7,
                       show_xticks = TRUE, ...) {
  order_cells <- match.arg(order_cells)
  ord_table <- NULL
  if (inherits(x, "deconv")) {
    scheme <- material_colours(x$mk)
    ord_table <- order(x$mk$cell_table)
    x <- x$subclass$output
  }
  if (is.null(scheme)) scheme <- hue_pal(h = c(0, 270))(ncol(x))
  
  if (length(order_col) > 1 || order_col != 0) {
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
    x <- x[ord, ]
  }
  
  if (!is.null(seriate)) {
    if (!requireNamespace("seriation", quietly = TRUE)) {
      stop("Package 'seriation' is not installed", call. = FALSE)
    }
    ord <- seriation::seriate(x, method = seriate, margin = 1)[[1]]
    x <- x[ord, ]
  }
  
  cell_ord <- order(colMeans(x))
  if (order_cells == "increase") {
    x <- x[, cell_ord]
  } else if (order_cells == "decrease") {
    x <- x[, rev(cell_ord)]
  } else if (!is.null(ord_table)) {
    x <- x[, ord_table]
  }
  
  strw <- max(strwidth(rownames(x), units = "inches", cex = cex.names),
              na.rm = TRUE)
  mar1 <- strw / par("csi") +1.5
  if (!show_xticks) mar1 <- 2
  op <- par(mar = c(mar1, 4, 1.5, 1.5))
  on.exit(par(op))
  barplot(t(x), las = 2, col = scheme,
          cex.names = cex.names, xaxt = ifelse(show_xticks, "s", "n"),
          tcl = -0.3, mgp = c(2.2, 0.5, 0), ...)
  # extend axis line
  yrange <- par("usr")[3:4]
  pos <- par("usr")[1]
  axis(2, at = c(0, yrange[2]), pos = pos, labels = FALSE, lwd.ticks = 0)
}

#' @rdname stack_plot
#' @importFrom ggplot2 ggplot geom_col aes scale_fill_manual xlab ylab
#'   theme_classic theme element_text guide_legend guides guide_axis
#'   scale_y_continuous expansion element_blank element_line
#' @importFrom rlang .data
#' @importFrom utils stack
#' @export

stack_ggplot <- function(x, percent = FALSE, order_col = 1, scheme = NULL,
                         order_cells = c("none", "increase", "decrease"),
                         seriate = NULL,
                         legend_ncol = NULL, legend_position = "bottom",
                         show_xticks = FALSE) {
  ord_table <- NULL
  if (inherits(x, "deconv")) {
    scheme <- material_colours(x$mk)
    ord_table <- order(x$mk$cell_table)
    x <- x$subclass$output
  }
  if (is.null(scheme)) {
    scheme <- hue_pal(h = c(0, 270))(ncol(x))
  }
  if (percent) {
    rs <- rowSums(x)
    x <- x / rs * 100
  }
  order_cells <- match.arg(order_cells)
  ord <- seq_len(nrow(x))
  if (length(order_col) > 1 || order_col != 0) {
    if (percent) {
      ord <- if (length(order_col) == 1) {
        order(x[, order_col])
      } else {
        order(rowMeans(x[, order_col, drop = FALSE]))
      }
    } else {
      ord <- order(rowSums(x))
    }
  }
  
  if (!is.null(seriate)) {
    if (!requireNamespace("seriation", quietly = TRUE)) {
      stop("Package 'seriation' is not installed", call. = FALSE)
    }
    ord <- seriation::seriate(x, method = seriate, margin = 1)[[1]]
  }
  
  ylab <- if (percent) "Cell proportion(%)" else "Relative cell amount"
  
  df <- stack(as.data.frame(x))
  df$id <- factor(rep(rownames(x), ncol(x)), levels = rownames(x)[ord])
  
  cell_ord <- order(colMeans(x))
  if (order_cells == "increase") {
    df$ind <- factor(df$ind, levels = colnames(x)[cell_ord])
  } else if (order_cells == "decrease") {
    df$ind <- factor(df$ind, levels = colnames(x)[rev(cell_ord)])
  } else if (!is.null(ord_table)) {
    df$ind <- factor(df$ind, levels = colnames(x)[ord_table])
  }
  
  p <- ggplot(df, aes(x = .data$id, y = .data$values, fill = .data$ind)) +
    geom_col(colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = scheme,
                      guide = guide_legend(title = "Cell type",
                                           ncol = legend_ncol,
                                           title.position = "top",
                                           position = legend_position)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    guides(x = guide_axis(angle = 90)) +
    xlab("") + ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.key.size = unit(0.8, 'lines'),
          legend.spacing.y = unit(0, 'lines'))
  if (!show_xticks) p <- p + theme(axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank())
  p
}


#' Cell subclass violin plot
#' 
#' Produces violin plots using ggplot2 showing amounts of cell subclasses in
#' deconvoluted bulk samples.
#' 
#' @param x matrix of deconvolution results with samples in rows and cell
#'   subclasses or groups in columns. If a 'deconv' class object is supplied the
#'   deconvolution values for the cell subclasses are extracted and plotted.
#' @param percent Logical whether to scale the matrix rows as percentage.
#' @param order_cols Character value specifying with cell types are ordered by
#'   mean abundance.
#' @returns A ggplot2 plotting object.
#' @importFrom ggplot2 geom_violin
#' @export
violin_plot <- function(x, percent = FALSE,
                        order_cols = c("none", "increase", "decrease")) {
  order_cols <- match.arg(order_cols)
  if (inherits(x, "deconv")) {
    x <- if (percent) x$subclass$percent else x$subclass$output
  }
  if (order_cols != "none") {
    o <- order(colMeans(x), decreasing = (order_cols == "decrease"))
    x <- x[, o]
  }
  ylab <- if (percent) "Cell proportion(%)" else "Relative cell amount"
  df <- stack(as.data.frame(x))
  
  ggplot(df, aes(x = .data$ind, y = .data$values, fill = .data$ind)) +
    geom_violin(scale = "width") +
    guides(x = guide_axis(angle = 60)) +
    xlab("") + ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          legend.position = "none")
}


material_colours <- function(mk) {
  subcl_tab <- table(mk$cell_table)
  ngroup <- length(subcl_tab)
  if (ngroup < 20 && requireNamespace("ggsci", quietly = TRUE)) {
    pal <- eval(formals(ggsci::pal_material)$palette)
    group_pal <- pal[round(19 / ngroup * seq(0, ngroup-1)) +1]
    cols <- lapply(seq_len(ngroup), function(i) {
      ggsci::pal_material(group_pal[i], subcl_tab[i])(subcl_tab[i])
    })
    return(unlist(cols))
  }
  hue_pal(h = c(0, 300))(length(mk$cell_table))
}
