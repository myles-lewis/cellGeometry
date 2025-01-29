
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
#' used with [deconvolute()]: `exp_signature = TRUE`, `convert_bulk = FALSE`,
#' `use_filter = FALSE` and `comp_amount = 1`.
#' 
#' @param object Either a 'cellMarkers' class object or a single cell count
#'   matrix with genes in rows and cells in columns.
#' @param samples An integer matrix with samples in rows and columns for each
#'   cell subclass in `object`. This can be generated using [generate_samples()].
#' @param subclass Vector of cell subclasses matching the columns in `object`.
#'   Only used if `object` is a single cell count matrix.
#' @param times Scaling factor to increase sampling of cells. Each cell in
#'   `samples` is randomly sampled this many times. Only used if `object` is a
#'   single cell count matrix.
#' @returns An integer count matrix with genes in rows and cell subclasses in
#'   columns. This can be used as `test` with the [deconvolute()] function.
#' @seealso [generate_samples()] [deconvolute()]
#' @export
simulate_bulk <- function(object, samples, subclass, times = 300) {
  if (inherits(object, "cellMarkers")) {
    genemean_counts <- 2^object$genemeans -1
    if (ncol(genemean_counts) != ncol(samples)) stop("incompatible number of columns")
    sim_pseudo <- genemean_counts %*% t(samples)
    mode(sim_pseudo) <- "integer"
    return(sim_pseudo)
  }
  if (!inherits(object, c("dgCMatrix", "matrix", "Seurat"))) object <- as.matrix(object)
  if (ncol(object) != length(subclass)) stop("incompatible dimensions")
  if (!is.factor(subclass)) subclass <- factor(subclass)
  if (any(!colnames(samples) %in% levels(subclass))) stop("incompatible subclasses")
  message("Creating sampling matrix")
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
  message("Matrix multiplication")
  sim_pseudo <- as.matrix(object %*% cmat)
  colnames(sim_pseudo) <- rownames(samples)
  if (max(sim_pseudo) <= .Machine$integer.max) mode(sim_pseudo) <- "integer"
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
#' @param force_intercept Logical whether to force intercept through 0.
#' @param show_identity Logical whether to show the identity line.
#' @param cols Optional vector of column indices to plot to show either a subset
#'   of columns or change the order in which columns are plotted.
#' @param title Title for page of plots.
#' @param cex.title Font size for title.
#' @param ... Optional arguments passed to `plot()`.
#' @returns No return value. Produces scatter plots using base graphics.
#' @importFrom graphics abline mtext
#' @importFrom stats lm runif
#' @export
plot_set <- function(obs, pred, mfrow = NULL,
                     force_intercept = FALSE,
                     show_identity = FALSE,
                     cols = NULL,
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
  col <- if (force_intercept) "red" else "blue"
  for (i in subclasses) {
    plot(obs[, i], pred[, i], cex = 0.8, pch = 16,
         xlab = i, ylab = "Predicted", ...)
    fit <- if (force_intercept) {
      lm(pred[, i] ~ obs[, i] + 0)
    } else lm(pred[, i] ~ obs[, i])
    abline(fit, col = col)
    if (show_identity) abline(0, 1, col = "grey50", lty = 2)
    rsq <- summary(fit)$r.squared |> format(digits = 3)
    mtext(bquote(R^2 == .(rsq)), cex = par("cex"), adj = 0.04)
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
