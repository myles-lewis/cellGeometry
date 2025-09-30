
#' Specificity plot
#' 
#' Scatter plot showing specificity of genes as markers for a particular cell
#' subclass. Optimal gene markers for that cell subclass are those genes which
#' are closest to or lie on the y axis, while also being of highest mean
#' expression.
#' 
#' For `type = 1`, coordinates are drawn as x = angle of vector in degrees, y =
#' mean gene expression of each gene in the subclass of interest. This version
#' is easier to use to identify additional gene markers. The plotly version
#' allows users to hover over points and identify which gene they belong to.
#' 
#' If `type = 2`, the coordinates are drawn as x = vector length * sin(angle)
#' and y = vector length * cos(angle), where vector length is the Euclidean
#' length of that gene in space where each cell subclass is a dimension. Angle
#' is the angle between the projected vector in space against perfection for
#' that cell subclass, i.e. the vector lying perfectly along the subclass
#' dimension with no deviation along other subclass dimensions, i.e. a gene
#' which is expressed solely in that subclass and has 0 expression in all other
#' subclasses. y is equal to the mean expression of each gene in the subclass of
#' interest. x represents the Euclidean distance of mean expression in all other
#' subclasses, i.e. overall non-specific gene expression in other subclasses.
#' Thus, the plot represents a rotation of all genes as vectors around the axis
#' of the subclass of interest onto the same plane so that the angle with the
#' subclass of interest is visualised between genes.
#' 
#' Colour is used to overlay the ranking of each gene across the subclasses,
#' showing for each gene where the subclass of interest is ranked compared to
#' the other subclasses. Best markers have the subclass of interest ranked 1st.
#' 
#' @param mk a 'cellMarkers' class object.
#' @param subclass character value specifying the subclass to be plotted.
#' @param group character value specifying cell group to be plotted. One of
#'   `subclass` or `group` must be specified.
#' @param type Numeric value, either 1 (the default) for a plot of angle on x
#'   axis and mean expression on y axis; or 2 for a plot projecting the vector
#'   angle into the same plain. See Details below.
#' @param use_filter logical, whether to use gene mean expression to which
#'   noise reduction filtering has been applied.
#' @param nrank number of ranks of subclasses to display.
#' @param nsubclass numeric value, number of top markers to label. By default
#'   this is obtained from `mk` for that subclass.
#' @param expfilter numeric value for the expression filter level below which
#'   genes are excluded from being markers. Defaults to the level used when
#'   `cellMarkers()` or `updateMarkers()` was called.
#' @param scheme Vector of colours for points.
#' @param add_labels character vector of additional genes to label
#' @param label_pos character value, either "left" or "right" specifying which
#'   side to add labels. Only for `type = 1` plots.
#' @param axis_extend numeric value, specifying how far to extend the x axis to
#'   the left as a proportion. Only invoked when `label_pos = "left"`.
#' @param nudge_x,nudge_y Label adjustments passed to `geom_label_repel()` or
#'   `geom_text_repel()`.
#' @param ... Optional arguments passed to `geom_label_repel()` or
#'   `geom_text_repel()` for `specificity_plot()` or `plot_ly()` for
#'   `specificity_plotly()`.
#' @returns ggplot2 or plotly scatter plot object.
#' @importFrom ggplot2 geom_point geom_vline scale_color_manual xlim ylim
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom grDevices adjustcolor
#' @export

specificity_plot <- function(mk, subclass = NULL,
                             group = NULL,
                             type = 1,
                             use_filter = FALSE,
                             nrank = 8,
                             nsubclass = NULL,
                             expfilter = NULL,
                             scheme = NULL,
                             add_labels = NULL,
                             label_pos = "right",
                             axis_extend = 0.4,
                             nudge_x = NULL, nudge_y = NULL,
                             ...) {
  if (!inherits(mk, "cellMarkers")) stop("not a 'cellMarkers' class object")
  if (is.null(subclass) & is.null(group))
    stop("Either subclass or group must be specified")
  
  if (!is.null(subclass)) {
    if (is.numeric(subclass)) subclass <- colnames(mk$genemeans)[subclass]
    if (!subclass %in% colnames(mk$genemeans))
      stop("subclass ", subclass, " not found")
    genemeans <- if (use_filter) mk$genemeans_filtered else mk$genemeans
    if (is.null(nsubclass)) nsubclass <- mk$nsubclass[subclass]
    if (is.null(nsubclass)) nsubclass <- 5
    labs <- rownames(mk$best_angle[[subclass]][1L:nsubclass, ])
    subc <- subclass
  } else {
    if (is.numeric(group)) group <- colnames(mk$groupmeans)[group]
    if (!group %in% colnames(mk$groupmeans))
      stop("group ", group, " not found")
    genemeans <- if (use_filter) mk$groupmeans_filtered else mk$groupmeans
    if (is.null(nsubclass)) nsubclass <- 5
    labs <- rownames(mk$group_angle[[group]][1L:nsubclass, ])
    subc <- group
  }
  
  vecLength <- sqrt(rowSums(genemeans^2))
  genemeans_scaled <- genemeans / vecLength
  genemeans_angle <- acos(genemeans_scaled)
  gene_rank <- apply(-genemeans, 1, rank)[subc, ]
  nrank <- pmin(ncol(genemeans), nrank)
  gene_rank[gene_rank > nrank] <- nrank
  gene_rank <- factor(floor(gene_rank), levels = 1:nrank)
  if (type == 1) {
    if (is.null(expfilter)) expfilter <- mk$opt$expfilter
    low <- genemeans[, subc] < expfilter & !rownames(genemeans) %in% labs
    gene_rank[low] <- nrank
    levels(gene_rank)[nrank] <- paste0(nrank, "+/low")
  } else {
    levels(gene_rank)[nrank] <- paste0(nrank, "+")
  }
  
  df <- data.frame(angle = genemeans_angle[, subc],
                   angle.deg = genemeans_angle[, subc] * 180/pi,
                   mean = genemeans[, subc],
                   rank = gene_rank)
  df$x <- vecLength * sin(df$angle)
  df$y <- vecLength * cos(df$angle)
  df <- df[vecLength != 0, ]
  labs <- unique(c(labs, add_labels))
  df$label <- ""
  df$label[match(labs, rownames(df))] <- labs
  df <- df[rev(order(df$rank)), ]
  
  if (is.null(scheme)) {
    scheme <- c(hue_pal(h = c(0, 270), c = 120)(nrank -1),
                adjustcolor("grey", 0.5))
    scheme[1] <- "red"
  }
  
  if (type == 2) {
    # use actual angle; radius is vecLength
    xlim <- xr <- range(df$x, na.rm = TRUE)
    yr <- range(df$y, na.rm = TRUE)
    if (label_pos == "left") {
      xlim[1] <- xlim[1] - diff(xr) * axis_extend
      if (is.null(nudge_x)) nudge_x <- -diff(xr) * 0.2
    } else {
      if (is.null(nudge_x)) nudge_x <- diff(xr) * 0.5
    }
    if (is.null(nudge_y)) nudge_y <- 0
    
    ggplot(df, aes(x = .data$x, y = .data$y, color = .data$rank,
                   label = .data$label)) +
      (if (label_pos == "left") geom_vline(xintercept = 0, lty = 2)) +
      geom_point(show.legend = TRUE) +
      scale_color_manual(values = scheme, drop = FALSE) +
      (if (label_pos == "left") {
        geom_label_repel(size = 3, color = "black",
                         nudge_x = nudge_x, nudge_y = nudge_y,
                         hjust = 1, label.size = NA, direction = "y", ...)
      } else {
        geom_text_repel(size = 3, color = "black",
                        nudge_x = nudge_x, nudge_y = nudge_y,
                        hjust = 0, direction = "y", na.rm = TRUE, ...)
      }) +
      xlim(xlim) + ylim(yr) +
      xlab("Non-specific gene expression") +
      ylab(paste(subc, "mean expression")) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(color = "black"))
  } else {
    # angle on x, mean exp on y
    xr <- range(df$angle.deg, na.rm = TRUE)
    yr <- range(df$mean, na.rm = TRUE)
    if (is.null(nudge_x)) nudge_x <- 0
    if (is.null(nudge_y)) nudge_y <- 0.1
    
    ggplot(df, aes(x = .data$angle.deg, y = .data$mean, color = .data$rank,
                   label = .data$label)) +
      geom_point(show.legend = TRUE) +
      scale_color_manual(values = scheme, drop = FALSE) +
      geom_text_repel(size = 3, color = "black",
                      nudge_x = nudge_x, nudge_y = nudge_y, ...) +
      xlab("Vector angle") +
      ylab(paste(subc, "mean expression")) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(color = "black"))
  }
}


#' @rdname specificity_plot
#' @export
specificity_plotly <- function(mk, subclass = NULL,
                               group = NULL,
                               type = 1,
                               use_filter = FALSE,
                               nrank = 8,
                               nsubclass = NULL,
                               expfilter = NULL,
                               scheme = NULL,
                               ...) {
  if (!inherits(mk, "cellMarkers")) stop("not a 'cellMarkers' class object")
  if (is.null(subclass) & is.null(group))
    stop("Either subclass or group must be specified")
  
  if (!is.null(subclass)) {
    if (is.numeric(subclass)) subclass <- colnames(mk$genemeans)[subclass]
    if (!subclass %in% colnames(mk$genemeans))
      stop("subclass ", subclass, " not found")
    genemeans <- if (use_filter) mk$genemeans_filtered else mk$genemeans
    if (is.null(nsubclass)) nsubclass <- mk$nsubclass[subclass]
    if (is.null(nsubclass)) nsubclass <- 5
    labs <- rownames(mk$best_angle[[subclass]][1L:nsubclass, ])
    subc <- subclass
  } else {
    if (is.numeric(group)) group <- colnames(mk$groupmeans)[group]
    if (!group %in% colnames(mk$groupmeans))
      stop("group ", group, " not found")
    genemeans <- if (use_filter) mk$groupmeans_filtered else mk$groupmeans
    if (is.null(nsubclass)) nsubclass <- 5
    labs <- rownames(mk$group_angle[[group]][1L:nsubclass, ])
    subc <- group
  }
  
  vecLength <- sqrt(rowSums(genemeans^2))
  genemeans_scaled <- genemeans / vecLength
  genemeans_angle <- acos(genemeans_scaled)
  gene_rank <- apply(-genemeans, 1, rank)[subc, ]
  nrank <- pmin(ncol(genemeans), nrank)
  gene_rank[gene_rank > nrank] <- nrank
  gene_rank <- factor(floor(gene_rank))
  if (type == 1) {
    if (is.null(expfilter)) expfilter <- mk$opt$expfilter
    low <- genemeans[, subc] < expfilter & !rownames(genemeans) %in% labs
    gene_rank[low] <- nrank
    levels(gene_rank)[nrank] <- paste0(nrank, "+/low")
  } else {
    levels(gene_rank)[nrank] <- paste0(nrank, "+")
  }
  gene_rank <- factor(gene_rank, levels = rev(levels(gene_rank)))
  
  df <- data.frame(angle = genemeans_angle[, subc],
                   angle.deg = genemeans_angle[, subc] * 180/pi,
                   mean = genemeans[, subc],
                   rank = gene_rank)
  df$x <- vecLength * sin(df$angle)
  df$y <- vecLength * cos(df$angle)
  df <- df[vecLength != 0, ]
  
  df$text <- paste(rownames(df), "<br>Angle: ", signif(df$angle.deg, 3),
                   "<br>Mean expr:", signif(df$mean, 3),
                   "<br>Rank:", df$rank)
  
  if (is.null(scheme)) {
    scheme <- c(hue_pal(h = c(0, 270), c = 120)(nrank -1),
                adjustcolor("grey", 0.5))
    scheme[1] <- "red"
  }
  scheme <- rev(scheme)
  
  if (!requireNamespace("plotly", quietly = TRUE)) 
    stop("Package 'plotly' is not installed", call. = FALSE)
  if (type == 2) {
    # use actual angle; radius is vecLength
    plotly::plot_ly(df, x = ~x, y = ~y, color = ~rank, colors = scheme,
            mode = "markers", type = "scattergl", hoverinfo = "text",
            marker = list(size = 6.5, line = list(width = 0.2, color = "white")),
            text = ~text, ...) |>
      plotly::layout(legend = list(traceorder = "reversed"),
             xaxis = list(title = "Non-specific gene expression"),
             yaxis = list(title = paste(subc, "mean expression")))
  } else {
    # angle on x, mean exp on y
    plotly::plot_ly(df, x = ~angle.deg, y = ~mean, color = ~rank,
            colors = scheme,
            mode = "markers", type = "scattergl", hoverinfo = "text",
            marker = list(size = 6.5, line = list(width = 0.2, color = "white")),
            text = ~text, ...) |>
      plotly::layout(legend = list(traceorder = "reversed"),
           xaxis = list(title = "Angle"),
           yaxis = list(title = paste(subc, "mean expression")))
  }
}
