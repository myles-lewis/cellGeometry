
#' Residuals plot
#' 
#' Plots residuals from a deconvolution result object against bulk gene
#' expression (on semi-log axis). Normal residuals, weighted residuals or
#' Studentized residuals can be visualised to check for heteroscedasticity and
#' genes with extreme errors.
#' 
#' @param fit 'deconv' class deconvolution object
#' @param test bulk gene expression matrix assumed to be in raw counts
#' @param type Specifies type of residuals to be plotted
#' @param show_outliers Logical whether to show any remaining outlying extreme
#'   genes in red
#' @param show_plot Logical whether to show plot using base graphics (used to
#'   allow return of dataframe of points without plotting)
#' @param nlabels Number of outlying gene labels to add. Labels are only added
#'   if outliers are detected. Although each gene has a point for the residual
#'   for every bulk sample, only the point with the maximum absolute residual is
#'   labelled.
#' @param ... Optional arguments passed to [plot()] for the base graphics
#'   version, or to `geom_text_repel()` for the ggplot2 version.
#' @returns Produces a scatter plot in base graphics. Returns invisibly a
#'   dataframe of the coordinates of the points. The ggplot version returns a
#'   ggplot2 plotting object.
#' @export
plot_residuals <- function(fit, test, type = c("raw", "student", "weight"),
                           show_outliers = TRUE, show_plot = TRUE,
                           ...) {
  if (!inherits(fit, "deconv")) stop("not a 'deconv' class object")
  type <- match.arg(type)
  geneset <- rownames(fit$subclass$residuals)
  ge <- test[geneset, ]
  res <- fit$subclass$residuals
  w <- fit$subclass$weights %||% 1
  res <- switch(type, student = rstudent(fit),
                weight = res * w,
                res)
  ylab <- switch(type, student = "Studentized residuals",
                 weight = "Weighted residuals",
                 "Raw residuals")
  col <- adjustcolor("black", 0.2)
  dat <- data.frame(expr = as.vector(ge), res = as.vector(res),
                    gene = rownames(res))
  if (show_outliers) {
    outlier_method <- fit$opt$outlier_method %||% "var.e"
    outlier_cutoff <- fit$opt$outlier_cutoff %||%
      switch(outlier_method, var.e = 4, cooks = 1, rstudent = 10)
    outlier_quantile <- fit$opt$outlier_quantile %||% 0.9
    count_space <- fit$opt$count_space %||% TRUE
    metric <- outlier_metric(fit$subclass, outlier_method, outlier_quantile,
                             count_space)
    outlier <- metric > outlier_cutoff
    nm <- names(metric[outlier])
    dat$outlier <- FALSE
    dat$outlier[dat$gene %in% nm] <- TRUE
    dat <- dat[order(dat$outlier), ]
    col <- rep(col, nrow(dat))
    col[dat$outlier] <- adjustcolor("red", 0.7)
  }
  if (show_plot) {
    nm <- fit$call$test
    .call <- match.call()
    if (nm != .call$test) {
      message("`", .call$test, "` does not match the bulk dataset `", nm,
              "` used in `", .call$fit, "`")
    }
    new.args <- list(...)
    args <- list(x = dat$expr, y = dat$res, log = "x",
                 xlim = c(pmax(1, min(ge, na.rm = TRUE)), max(ge, na.rm = TRUE)),
                 cex = 0.7, col = col, pch = 16, bty = "l",
                 xlab = "", ylab = ylab)
    if (length(new.args)) args[names(new.args)] <- new.args     
    do.call("plot", args) |>
      suppressWarnings()
    if (is.null(new.args$xlab)) {
      mtext("Bulk gene expression", 1, line = 2.2, cex = par("cex.lab") * par("cex"))
    }
    abline(0, 0, col = "royalblue")
  }
  invisible(dat)
}


#' @rdname plot_residuals
#' @importFrom ggplot2 scale_x_log10 scale_colour_manual geom_hline
#' @export
ggplot_residuals <- function(fit, test, type = c("reg", "student", "weight"),
                             show_outliers = TRUE, nlabels = 0, ...) {
  nm <- fit$call$test
  .call <- match.call()
  if (nm != .call$test) {
    message("`", .call$test, "` does not match the bulk dataset `", nm,
            "` used in `", .call$fit, "`")
  }
  dat <- plot_residuals(fit, test, type, show_outliers, show_plot = FALSE)
  type <- match.arg(type)
  dat$expr[dat$expr < 1] <- NA
  ylab <- switch(type, student = "Studentized residuals",
                 weight = "Weighted residuals",
                 "Raw residuals")
  dat$label <- ""
  outl <- unique(dat$gene[dat$outlier])
  nlabels <- min(c(nlabels, length(outl)))
  if (nlabels > 0) {
    m <- vapply(outl, function(i) {
      abs(mean(dat$res[dat$gene == i]))
    }, numeric(1))
    outl <- outl[order(m, decreasing = TRUE)]
    labs <- outl[seq_len(nlabels)]
    for (i in labs) {
      ind <- which(dat$gene == i)
      w <- which.max(dat$res[ind])
      dat$label[ind[w]] <- dat$gene[ind[w]]
    }
  }
  
  ggplot(dat, aes(x = .data$expr, y = .data$res, label = .data$label)) +
    geom_point(aes(colour = .data$outlier), na.rm = TRUE) +
    geom_hline(yintercept = 0, color = "royalblue") +
    scale_x_log10(guide = "axis_logticks") +
    scale_colour_manual(values = c(adjustcolor("black", 0.2),
                                   adjustcolor("red", 0.7))) +
    (if (nlabels > 0) {
      geom_text_repel(size = 3, color = "black", na.rm = TRUE, ...)
    }) +
    xlab("Bulk gene expression") + ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none")
}
