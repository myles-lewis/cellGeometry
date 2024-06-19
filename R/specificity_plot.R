
#' Specificity plot
#' 
#' Scatter plot showing specificity of genes as markers for a particular cell
#' subclass. Optimal gene markers for that cell subclass are those genes which
#' are closest to or lie on the y axis, while also being of highest mean
#' expression.
#' 
#' The coordinates are drawn as x = vector length * sin(angle) and y = vector
#' length * cos(angle), where vector length is the Euclidean length of that gene
#' in space where each cell subclass is a dimension. Angle is the angle between
#' the projected vector in space against perfection for that cell
#' subclass, i.e. the vector lying perfectly along the subclass dimension with
#' no deviation along other subclass dimensions, i.e. a gene which is expressed
#' solely in that subclass and has 0 expression in all other subclasses. y is
#' equal to the mean expression of each gene in the subclass of interest. x
#' represents the Euclidean distance of mean expression in all other subclasses.
#' 
#' Colour is used to overlay the ranking of each gene across the subclasses,
#' showing for each gene where the subclass of interest is ranked compared to
#' the other subclasses. Best markers have the subclass of interest ranked 1st.
#' 
#' @param mk a 'cellMarkers' class object.
#' @param subclass character value specifying the subclass to be plotted
#' @param use_filter logical, whether to use gene mean expressions to which
#'   noise reduction filtering has been applied.
#' @param nrank number of ranks of subclasses to display.
#' @param label_pos character value, either "left" or "right" specifying which
#'   side to add labels.
#' @param nsubclass number of top markers to label
#' @param axis_extend numeric value, specifying how far to extend the x axis to
#'   the left as a proportion. Only invoked when `label_pos = "left"`.
#' @param scheme Vector of colours for points.
#' @returns ggplot2 scatter plot object.
#' @importFrom ggplot2 geom_point geom_vline scale_color_manual xlim ylim
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom grDevices adjustcolor
#' @export

specificity_plot <- function(mk, subclass,
                             use_filter = FALSE,
                             nrank = 8,
                             label_pos = "right",
                             nsubclass = 5,
                             axis_extend = 0.4,
                             scheme = NULL) {
  if (!inherits(mk, "cellMarkers")) stop("not a 'cellMarkers' class object")
  if (is.numeric(subclass)) subclass <- colnames(mk$genemeans)[subclass]
  if (!subclass %in% colnames(mk$genemeans))
    stop("subclass ", subclass, " not found")
  
  genemeans <- if (use_filter) mk$genemeans_filtered else mk$genemeans
  vecLength <- sqrt(rowSums(genemeans^2))
  genemeans_scaled <- genemeans / vecLength
  genemeans_angle <- acos(genemeans_scaled)
  gene_rank <- apply(-genemeans, 1, rank)[subclass, ]
  gene_rank[gene_rank > nrank] <- nrank
  gene_rank <- factor(floor(gene_rank))
  levels(gene_rank)[nrank] <- paste0(nrank, "+")
  
  df <- data.frame(angle = genemeans_angle[, subclass],
                   angle.deg = genemeans_angle[, subclass] * 180/pi,
                   mean = genemeans[, subclass],
                   rank = gene_rank)
  df$x <- vecLength * sin(df$angle)
  df$y <- vecLength * cos(df$angle)
  df <- df[vecLength != 0, ]
  labs <- rownames(mk$best_angle[[subclass]][1:nsubclass, ])
  df$label <- NA
  df$label[match(labs, rownames(df))] <- labs
  df <- df[rev(order(df$rank)), ]
  
  if (is.null(scheme)) {
    scheme <- c(hue_pal(h = c(0, 270), c = 120)(nrank -1),
                adjustcolor("grey", 0.5))
    scheme[1] <- "red"
  }
  
  xlim <- xr <- range(df$x, na.rm = TRUE)
  yr <- range(df$y, na.rm = TRUE)
  if (label_pos == "left") xlim[1] <- xlim[1] - diff(xr) * axis_extend
  
  ## use actual angle; radius is vecLength
  ggplot(df, aes(x = .data$x, y = .data$y, color = .data$rank,
                 label = .data$label)) +
    (if (label_pos == "left") geom_vline(xintercept = 0, lty = 2)) +
    geom_point() +
    scale_color_manual(values = scheme) +
    (if (label_pos == "left") {
      geom_label_repel(size = 3, color = "black", nudge_x = -diff(xr) * 0.2,
                       hjust = 1, label.size = NA, direction = "y", na.rm = TRUE)
    } else {
      geom_text_repel(size = 3, color = "black", nudge_x = diff(xr) * 0.5,
                      hjust = 0, direction = "y", na.rm = TRUE)
    }) +
    xlim(xlim) + ylim(yr) +
    xlab("Vector length * sin(angle)") +
    ylab(paste(subclass, "mean expression")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))
}
