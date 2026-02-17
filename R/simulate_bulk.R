
#' Generate random cell number samples
#' 
#' Used for simulating pseudo-bulk RNA-Seq from a 'cellMarkers' object. Cell
#' counts are randomly sampled from the uniform distribution, using the original
#' subclass contingency table as a limit on the maximum number of cells in each
#' subclass.
#' 
#' @param object A 'cellMarkers' class object
#' @param n Integer value for the number of samples to generate
#' @param equal_sample Logical whether to sample subclasses equally or generate
#'   samples with proportions of cells in keeping with the original subtotal of
#'   cells in the main scRNA-Seq data.
#' @param method Either "unif" or "dirichlet" to specify whether cell numbers
#'   are drawn from uniform distribution or dirichlet distribution.
#' @param alpha Shape parameter for `gtools::rdirichlet()`. Automatically
#'   expanded to be a vector whose length is the number of subclasses.
#' @param zero_fraction Numeric from 0 to 1 specifying proportion of cell
#'   subclasses to randomly set to zero in each sample.
#' @returns An integer matrix with `n` rows, with columns for each cell
#'   subclasses in `object`, representing cell counts for each cell subclass.
#'   Designed to be passed to [simulate_bulk()].
#' @details
#' Leaving `equal_sample = TRUE` is better for tuning deconvolution parameters.
#' 
#' @seealso [simulate_bulk()]
#' @importFrom gtools rdirichlet
#' @export
generate_samples <- function(object, n, equal_sample = TRUE,
                             method = c("unif", "dirichlet"),
                             alpha = 1.5, zero_fraction = 0) {
  lim <- object$subclass_table
  nc <- length(lim)
  method <- match.arg(method)
  if (method == "unif") {
    sim_counts <- matrix(runif(n * nc), ncol = nc,
                         dimnames = list(paste0("S", c(1:n)), names(lim)))
    if (equal_sample) {
      fac <- sum(lim) / nc
      sim_counts <- sim_counts * fac
    } else sim_counts <- t(t(sim_counts) * as.vector(lim))
  } else {
    sim_counts <- rdirichlet(n, rep_len(alpha, nc)) * sum(lim)
    dimnames(sim_counts) <- list(paste0("S", c(1:n)), names(lim))
  }
  mode(sim_counts) <- "integer"
  if (zero_fraction > 0) {
    zc <- round(zero_fraction * nc)
    wr <- rep(seq_len(n), each = zc)
    wc <- sample(nc, n * zc, replace = TRUE)
    w <- matrix(c(wr, wc), ncol = 2)
    sim_counts[w] <- 0L
  }
  sim_counts
}


#' Simulate pseudo-bulk RNA-Seq
#' 
#' Simulates pseudo-bulk RNA-Seq dataset using two modes. The first mode uses a
#' 'cellMarkers' class object and a matrix of counts for the numbers of cells of
#' each cell subclass. This method converts the log2 gene means back for
#' each cell subclass back to count scale and then calculates pseudo-bulk count
#' values based on the cell amounts specified in `samples`. In the 2nd mode, a
#' single-cell RNA-Seq dataset is required, such as a matrix used as input to
#' [cellMarkers()]. Cells from the relevant subclass are sampled from the
#' single-cell matrix in the appropriate amounts based on `samples`, except that
#' sampling is scaled up by the factor `times`.
#'
#' The first method can give perfect deconvolution if the following settings are
#' used with [deconvolute()]: `count_space = TRUE`, `convert_bulk = FALSE`,
#' `use_filter = FALSE` and `comp_amount = 1`.
#' 
#' @param object Either a 'cellMarkers' class object, or a single cell count
#'   matrix with genes in rows and cells in columns, with rownames representing
#'   gene IDs/symbols. The matrix can be a sparse matrix or DelayedMatrix.
#' @param samples An integer matrix of cell counts with samples in rows and
#'   columns for each cell subclass in `object`. This can be generated using
#'   [generate_samples()].
#' @param subclass Vector of cell subclasses matching the columns in `object`.
#'   Only used if `object` is a single cell count matrix.
#' @param times Scaling factor to increase sampling of cells. Cell counts in
#'   `samples` are scaled up by being multiplied by this number. Only used if
#'   `object` is a single cell count matrix.
#' @param method Either "dirichlet" or "unif" to specify whether cells are
#'   sampled based on the Dirichlet distribution with K = number of cells in
#'   each subclass, or sampled uniformly. When cells are oversampled uniformly,
#'   in the limit the summed gene expression tends to the arithmetic mean of the
#'   subclass x sample frequency. Dirichlet sampling provides proper randomness
#'   with sampling.
#' @param alpha Shape parameter for Dirichlet sampling.
#' @returns An integer count matrix with genes in rows and cell subclasses in
#'   columns. This can be used as `test` with the [deconvolute()] function.
#' @seealso [generate_samples()] [deconvolute()] [add_noise()]
#' @export
simulate_bulk <- function(object, samples, subclass, times = 1,
                          method = c("dirichlet", "unif"), alpha = 1) {
  if (inherits(object, "cellMarkers")) {
    genemean_counts <- 2^object$genemeans -1
    if (ncol(genemean_counts) != ncol(samples)) stop("incompatible number of columns")
    sim_pseudo <- genemean_counts %*% t(samples)
    mode(sim_pseudo) <- "integer"
    return(sim_pseudo)
  }
  
  # sample from count matrix
  if (!inherits(object, c("dgCMatrix", "matrix", "Seurat", "DelayedMatrix"))) {
    object <- as.matrix(object)
  }
  if (ncol(object) != length(subclass)) stop("incompatible dimensions")
  if (!is.factor(subclass)) subclass <- factor(subclass)
  if (any(!colnames(samples) %in% levels(subclass))) stop("incompatible subclasses")
  method <- match.arg(method)
  
  start <- Sys.time()
  message("Creating sampling matrix", appendLF = FALSE)
  samples <- samples * times
  subclass_lev <- levels(subclass)
  subclass_lev <- subclass_lev[subclass_lev %in% colnames(samples)]
  cmat <- vapply(seq_len(nrow(samples)), function(j) {
    s <- unlist(lapply(subclass_lev, function(i) {
      w <- which(subclass == i)
      n <- samples[j, i]
      if (method == "unif" || n < length(w) * 2) {
        sample(w, n, replace = TRUE)
      } else {
        rd <- as.vector(rdirichlet(1, rep_len(alpha, length(w))))
        sample(w, n, replace = TRUE, prob = rd)
      }
    }))
    tabulate(s, nbins = length(subclass))
  }, numeric(length(subclass)))
  cat_timer(start)
  
  start <- Sys.time()
  message("Matrix multiplication", appendLF = FALSE)
  sim_pseudo <- as.matrix(object %*% cmat)
  colnames(sim_pseudo) <- rownames(samples)
  if (max(sim_pseudo) <= .Machine$integer.max) mode(sim_pseudo) <- "integer"
  cat_timer(start)
  
  sim_pseudo
}

cat_timer <- function(start) {
  message(" (", format(Sys.time() - start, digits = 3), ")")
}

#' Scatter plots to compare deconvoluted subclasses
#' 
#' Produces a set of scatter plots using base graphics to compare actual cell
#' counts against deconvoluted cell counts from bulk (or pseudo-bulk) RNA-Seq
#' for each cell subclass. Mainly for use if ground truth is available, e.g. for
#' simulated pseudo-bulk RNA-Seq data.
#' 
#' @param obs Observed matrix of cell amounts with subclasses in columns and
#'   samples in rows.
#' @param pred Predicted (deconvoluted) matrix of cell amounts with rows and
#'   columns matching `obs`.
#' @param mfrow Optional vector of length 2 for organising plot layout. See
#'   `par()`.
#' @param show_zero Logical whether to force plot to include the origin.
#' @param show_identity Logical whether to show the identity line.
#' @param cols Optional vector of column indices to plot to show either a subset
#'   of columns or change the order in which columns are plotted. `NA` skips a
#'   plot space to introduce a gap between plots.
#' @param colour Colour for the regression lines.
#' @param title Title for page of plots.
#' @param cex.title Font size for title.
#' @param ... Optional arguments passed to `plot()`.
#' @returns No return value. Produces scatter plots using base graphics.
#' @importFrom graphics abline mtext plot.new
#' @importFrom stats lm runif rnorm
#' @export
plot_set <- function(obs, pred, mfrow = NULL,
                     show_zero = FALSE,
                     show_identity = FALSE,
                     cols = NULL,
                     colour = "blue",
                     title = "", cex.title = 1, ...) {
  if (!identical(dim(obs), dim(pred))) stop("incompatible dimensions")
  if (anyNA(pred)) {
    message("`pred` contains NA")
    pred[is.na(pred)] <- 0
  }
  if (is.null(cols)) cols <- TRUE
  subclasses <- colnames(obs)[cols]
  nr1 <- ceiling(sqrt(length(subclasses)))
  nr2 <- ceiling(length(subclasses) / nr1)
  if (is.null(mfrow)) mfrow <- c(nr1, nr2)
  oma <- par("oma")
  if (title != "" & oma[3] < 1.5) oma[3] <- 1.5
  op <- par(bty = "l", mgp = c(2.2, 0.6, 0), tcl = -0.3, oma = oma,
            mar = c(3.7, 3.7, 1.5, 1.1), mfrow = mfrow)
  on.exit(par(op))
  scheme <- rev(hue_pal(h = c(0, 270), c = 120)(11))
  xlim <- ylim <- NULL
  new.args <- list(...)
  
  for (i in subclasses) {
    if (is.na(i)) {plot.new(); next}
    if (show_zero) {
      xr <- range(obs[, i], na.rm = TRUE)
      xlim <- c(min(xr[1], 0), xr[2])
      yr <- range(pred[, i], na.rm = TRUE)
      ylim <- c(min(yr[1], 0), yr[2])
    }
    args <- list(x = obs[, i], y = pred[, i], cex = 0.8, pch = 16, las = 1,
                 xlab = i, ylab = "Predicted", xlim = xlim, ylim = ylim)
    if (length(new.args)) args[names(new.args)] <- new.args
    do.call(plot, args)
    fit <- lm(pred[, i] ~ obs[, i])
    rsq <- summary(fit)$r.squared
    col <- if (colour == "rainbow") scheme[ceiling(rsq*10) +1] else colour
    abline(fit, col = col, lwd = 1.5)
    if (show_identity) abline(0, 1, col = "grey50", lty = 2)
    mtext(bquote(R^2 == .(format(rsq, digits = 3))), cex = par("cex"), adj = 0.04)
  }
  mtext(title, outer = TRUE, cex = cex.title * par("cex"), adj = 0.05, line = 0)
}


#' Calculate R-squared and metrics on deconvoluted cell subclasses
#' 
#' Calculates Pearson r-squared, R-squared and RMSE comparing subclasses in each
#' column of `obs` with matching columns in deconvoluted `pred`. Samples are in
#' rows. For use if ground truth is available, e.g. simulated pseudo-bulk
#' RNA-Seq data.
#' 
#' Pearson r-squared ranges from 0 to 1. R-squared, calculated as 1 - rss/tss,
#' ranges from -Inf to 1.
#' 
#' @param obs Observed matrix of cell amounts with subclasses in columns and
#'   samples in rows.
#' @param pred Predicted (deconvoluted) matrix of cell amounts with rows and
#'   columns matching `obs`.
#' @returns Matrix containing Pearson r-squared, R-squared and RMSE values.
#' @importFrom stats cor
#' @export
metric_set <- function(obs, pred) {
  if (!identical(dim(obs), dim(pred))) stop("incompatible dimensions")
  if (anyNA(pred)) {
    message("`pred` contains NA")
    pred[is.na(pred)] <- 0
  }
  
  out <- t(vapply(colnames(obs), function(i) {
    r1 <- Rsq(obs[, i], pred[, i])
    r2 <- rmse(obs[, i], pred[, i])
    c(r1, r2)
  }, numeric(2)))
  
  cors <- diag(cor(obs, pred))^2 |> suppressWarnings()
  out <- cbind(cors, out)
  colnames(out) <- c("pearson.rsq", "Rsq", "RMSE")
  out
}

rmse <- function(obs, pred) {
  sqrt(mean((pred - obs)^2))
}

Rsq <- function(obs, pred) {
  rss <- sum((pred - obs)^2)
  tss <- sum((obs - mean(obs))^2)
  1 - rss/tss
}


#' Scatter plot to compare deconvoluted subclasses
#' 
#' Produces a single scatter plot using base graphics to compare actual cell
#' counts against deconvoluted cell counts from bulk (or pseudo-bulk) RNA-Seq
#' Cell subclasses are shown in different colours. Designed for use if ground
#' truth is available, e.g. for simulated pseudo-bulk RNA-Seq data.
#' 
#' @param obs Observed matrix of cell amounts with subclasses in columns and
#'   samples in rows.
#' @param pred Predicted (deconvoluted) matrix of cell amounts with rows and
#'   columns matching `obs`.
#' @param mk Optional matching cellMarkers object. This is used for its
#'   `cell_table` element to try to colour subclasses by group.
#' @param scheme Vector of colours, one for each cell subclass.
#' @param ellipse Either a single number for the number of ellipses to plot, in
#'   which case ellipses are shown for cell subclasses with the lowest R^2; or a
#'   character vector of cell subclasses to be outlined with ellipses. Requires
#'   the ggforce package to be installed.
#' @returns A ggplot2 scatter plot. An overall R^2 (coefficient of
#'   determination) comparing all observed and predicted results is shown.
#' @importFrom ggplot2 geom_abline alpha margin
#' @export
plot_pred <- function(obs, pred, mk = NULL, scheme = NULL, ellipse = NULL) {
  if (!identical(dim(obs), dim(pred))) stop("incompatible dimensions")
  if (anyNA(pred)) {
    message("`pred` contains NA")
    pred[is.na(pred)] <- 0
  }
  subclasses <- colnames(obs)
  if (is.null(scheme)) {
    if (is.null(mk)) {
      scheme <- hue_pal(h = c(0, 300))(ncol(obs))
    } else {
      if (!inherits(mk, "cellMarkers")) stop("mk is not a cellMarkers object")
      scheme <- material_colours(mk)
    }
  }
  dat <- data.frame(obs = as.vector(obs), pred = as.vector(pred),
                    subclass = factor(rep(subclasses, each = nrow(obs))))
  rsq <- format(Rsq(obs, pred), digits = 3)
  title <- bquote(R^2 ~"="~ .(rsq))
  
  p <- ggplot(dat, aes(x = .data$obs, y = .data$pred, color = .data$subclass,
                       fill = .data$subclass)) +
    geom_point(size = 1, alpha = 0.8) +
    scale_colour_manual(values = scheme) +
    scale_fill_manual(values = scheme) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Observed") + ylab("Predicted") +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          plot.title = element_text(size = 10),
          legend.position = "none")
  
  if (!is.null(ellipse)) {
    if (is.character(ellipse)) {
      subdat <- dat[dat$subclass %in% ellipse, ]
    } else {
      mset <- metric_set(obs, pred)
      o <- order(mset[, "Rsq"])[seq_len(ellipse)]
      subdat <- dat[dat$subclass %in% colnames(obs)[o], ]
    }
    if (!requireNamespace("ggforce", quietly = TRUE))
      stop("'ggforce' package is not installed", call. = FALSE)
    p <- p +
      ggforce::geom_mark_ellipse(data = subdat,
                                 aes(x = .data$obs, y = .data$pred,
                                     color = .data$subclass, fill = .data$subclass,
                                     label = .data$subclass),
                                 label.fontface = "plain",
                                 label.margin = margin(1, 1, 1, 1, "mm"),
                                 label.fontsize = 10,
                                 label.fill = alpha("white", 0.6),
                                 label.buffer = unit(5, "mm"),
                                 con.cap = 0,
                                 con.type = "straight",
                                 con.border = "none",
                                 expand = 0, linewidth = unit(0.3, "pt"))
  }
  p
}
