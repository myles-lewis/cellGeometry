
#' Spillover heatmap
#' 
#' Produces a heatmap from a 'cellMarkers' or 'deconv' class object showing
#' estimated amount of spillover between cell subclasses. The amount that each
#' cell subclass's overall vector spillovers (projects) into other cell
#' subclasses' vectors is shown in each row. Thus the column gives an estimate
#' of how much the most influential (specific) genes for a cell subclass are
#' expressed in other cells.
#' 
#' @param x Either a 'cellMarkers' or 'deconv' class object or a spillover
#'   matrix.
#' @param text Logical whether to show values of cells where spillover >
#'   `cutoff`. By default only shown for smaller matrices.
#' @param cutoff Threshold for showing values.
#' @param fontsize Numeric value for font size for cell values when
#'   `text = TRUE`.
#' @param subset Character vector of groups to be subsetted.
#' @param ... Optional arguments passed to [ComplexHeatmap::Heatmap()].
#' @returns No return value. Draws a heatmap using ComplexHeatmap.
#' @importFrom grid grid.text unit
#' @importFrom ComplexHeatmap draw
#' @importFrom circlize colorRamp2
#' @export

spillover_heatmap <- function(x,
                              text = NULL,
                              cutoff = 0.5,
                              fontsize = 8,
                              subset = NULL,
                              ...) {
  cell_table <- NULL
  if (inherits(x, "deconv")) {
    cell_table <- x$mk$cell_table
    x <- x$subclass$spillover
  }
  if (inherits(x, 'cellMarkers')) {
    cell_table <- x$cell_table
    m_itself <- x$spillover
  } else m_itself <- x
  
  if (!is.null(subset)) {
    s <- which(cell_table %in% subset)
    if (length(s) == 0) stop("no such subgroup")
    if (length(s) == 1) stop("subset too small")
    return(spillover_heatmap(x = m_itself[s, s], text = text, cutoff = cutoff,
                             fontsize = fontsize, ...))
  }
  
  mx <- max(m_itself)
  col <- if (mx > 1) {
    colorRamp2(c(0, 0.5, 0.8, 1, mx),
               c("#F4FAFF", "steelblue2", "purple", "red", "red3"))
  } else {
    colorRamp2(c(0, 0.5, 0.8, 1),
               c("#F4FAFF", "steelblue2", "purple", "red"))
  }
  
  if (is.null(text)) text <- dim(m_itself)[1] <= 50
  layer_fun <- NULL
  if (text) {
    layer_fun <- function(j, i, x, y, width, height, fill) {
      v <- pindex(m_itself, i, j)
      l <- v > cutoff & i != j
      if (any(l)) {
        grid.text(sprintf("%.1f", v[l]), x[l], y[l],
                  gp = gpar(fontsize = fontsize))
      }
    }
  }
  dots <- list(...)
  args <- list(m_itself, col = col,
               cluster_rows = FALSE, row_split = cell_table,
               cluster_row_slices = FALSE, row_title = NULL,
               cluster_columns = FALSE, column_split = cell_table,
               cluster_column_slices = FALSE, column_title = NULL,
               column_names_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               layer_fun = layer_fun,
               heatmap_legend_param = list(title = "spillover",
                                           legend_width = unit(6, "cm"),
                                           direction = "horizontal"))
  if (length(dots)) args[names(dots)] <- dots
  hm <- do.call(Heatmap, args) |> suppressMessages()
  draw(hm, heatmap_legend_side = "top")
}
