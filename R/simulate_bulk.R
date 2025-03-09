
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
#' @returns An integer matrix with `n` rows, with columns for each cell
#'   subclasses in `object`, representing cell counts for each cell subclass.
#'   Designed to be passed to [simulate_bulk()].
#' @details
#' Leaving `equal_sample = TRUE` is better for tuning deconvolution parameters.
#' 
#' @seealso [simulate_bulk()]
#' @export
generate_samples <- function(object, n, equal_sample = TRUE) {
  lim <- object$subclass_table
  nc <- length(lim)
  sim_counts <- matrix(runif(n * nc), ncol = nc,
                       dimnames = list(paste0("S", c(1:n)), names(lim)))
  if (equal_sample) {
    fac <- sum(lim) / nc
    sim_counts <- sim_counts * fac
  } else sim_counts <- t(t(sim_counts) * as.vector(lim))
  mode(sim_counts) <- "integer"
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
#' @param samples An integer matrix with samples in rows and columns for each
#'   cell subclass in `object`. This can be generated using [generate_samples()].
#' @param subclass Vector of cell subclasses matching the columns in `object`.
#'   Only used if `object` is a single cell count matrix.
#' @param times Scaling factor to increase sampling of cells. Each cell in
#'   `samples` is randomly sampled this many times. Only used if `object` is a
#'   single cell count matrix.
#' @returns An integer count matrix with genes in rows and cell subclasses in
#'   columns. This can be used as `test` with the [deconvolute()] function.
#' @seealso [generate_samples()] [deconvolute()] [add_noise()]
#' @export
simulate_bulk <- function(object, samples, subclass, times = 30,
                          add_noise = FALSE, sd1 = 100, sd2 = 0) {
  if (inherits(object, "cellMarkers")) {
    genemean_counts <- 2^object$genemeans -1
    if (ncol(genemean_counts) != ncol(samples)) stop("incompatible number of columns")
    sim_pseudo <- genemean_counts %*% t(samples)
    mode(sim_pseudo) <- "integer"
    return(sim_pseudo)
  }
  # sample from count matrix
  start <- Sys.time()
  if (!inherits(object, c("dgCMatrix", "matrix", "Seurat", "DelayedMatrix"))) {
    object <- as.matrix(object)
  }
  if (ncol(object) != length(subclass)) stop("incompatible dimensions")
  if (!is.factor(subclass)) subclass <- factor(subclass)
  if (any(!colnames(samples) %in% levels(subclass))) stop("incompatible subclasses")
  message("Creating sampling matrix", appendLF = FALSE)
  samples <- samples * times
  subclass_lev <- levels(subclass)
  subclass_lev <- subclass_lev[subclass_lev %in% colnames(samples)]
  cmat <- vapply(seq_len(nrow(samples)), function(j) {
    s <- unlist(lapply(subclass_lev, function(i) {
      w <- which(subclass == i)
      sample(w, samples[j, i], replace = TRUE)
    }))
    tabulate(s, nbins = length(subclass))
  }, numeric(length(subclass)))
  message(" (", format(Sys.time() - start, digits = 3), ")")
  start <- Sys.time()
  message("Matrix multiplication", appendLF = FALSE)
  sim_pseudo <- as.matrix(object %*% cmat)
  colnames(sim_pseudo) <- rownames(samples)
  if (max(sim_pseudo) <= .Machine$integer.max) mode(sim_pseudo) <- "integer"
  message(" (", format(Sys.time() - start, digits = 3), ")")
  
  sim_pseudo
}


#' Add noise to simulated count data
#' 
#' Gaussian noise can be added to the simulated count matrix in 2 ways which can
#' be combined. `sd1` controls addition of simple Gaussian noise to counts using
#' `rnorm` with sd specified by `sd1`. If `sd2` is >0, counts are converted
#' using log2+1 and Gaussian noise added, followed by conversion back to count
#' scale. Negative values are converted to 0. `sd1` affects low expressed genes
#' and hardly affects high expressed genes. `sd2` affects all genes irrespective
#' of expression level.
#' 
#' @param sim_pseudo An integer count matrix with genes in rows and cell
#'   subclasses typically generated by [simulate_bulk()].
#' @param sd1 Standard deviation of noise added to counts.
#' @param sd2 Standard deviation of noise added to log2(counts+1).
#' @returns An integer count matrix with genes in rows and cell subclasses in
#'   columns.
#' @export
add_noise <- function(sim_pseudo, sd1 = 100, sd2 = 0.05) {
  if (sd2 > 0) {
    # Gaussian noise on log scale
    rn <- rnorm(prod(dim(sim_pseudo)), sd = sd2)
    rmat <- matrix(rn, nrow = nrow(sim_pseudo))
    log_sim <- log2(sim_pseudo +1)
    log_sim <- log_sim + rmat
    sim_pseudo <- 2^log_sim -1
  }
  if (sd1 > 0) {
    # simple Gaussian noise
    rn <- rnorm(prod(dim(sim_pseudo)), sd = sd1)
    rmat <- matrix(round(rn), nrow = nrow(sim_pseudo))
    sim_pseudo <- sim_pseudo + rmat
  }
  
  sim_pseudo[sim_pseudo < 0] <- 0
  if (max(sim_pseudo) <= .Machine$integer.max) {
    mode(sim_pseudo) <- "integer"
  } else sim_pseudo <- round(sim_pseudo)
  sim_pseudo
}


#' Scatter plots to compare deconvoluted subclasses
#' 
#' Produces scatter plots using base graphics to compare actual cell counts
#' against deconvoluted cell counts from bulk (or pseudo-bulk) RNA-Seq. Mainly
#' for use if ground truth is available, e.g. for simulated pseudo-bulk RNA-Seq
#' data.
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
