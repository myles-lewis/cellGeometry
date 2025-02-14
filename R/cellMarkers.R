
#' Identify cell markers
#' 
#' Uses geometric method based on vector dot product to identify genes which are
#' the best markers for individual cell types.
#' 
#' If `verbose = TRUE`, the function will display an estimate of the required
#' memory. But importantly this estimate is only a guide. It is provided to help
#' users choose the optimal number of cores during parallelisation. Real memory
#' usage might well be more, theoretically up to double this amount, due to R's
#' use of copy-on-modify.
#' 
#' @param scdata Single-cell data matrix with genes in rows and cells in
#'   columns. Can be sparse matrix or DelayedMatrix. Must have rownames 
#'   representing gene IDs or gene symbols.
#' @param bulkdata Optional data matrix containing bulk RNA-Seq data with genes
#'   in rows. This matrix is only used for its rownames (gene IDs), to ensure
#'   that cell markers are selected from genes in the bulk dataset.
#' @param subclass Vector of cell subclasses matching the columns in `scdata`
#' @param cellgroup Optional grouping vector of major cell types matching the
#'   columns in `scdata`. `subclass` is assumed to contain subclasses which are
#'   subsets within `cellgroup` overarching classes.
#' @param nsubclass Number of genes to select for each single cell subclass.
#'   Either a single number or a vector with the number of genes for each
#'   subclass.
#' @param ngroup Number of genes to select for each cell group. Either a single
#'   number or a vector with the number of genes for each group.
#' @param expfilter Genes whose maximum mean expression on log2 scale per cell
#'   type are below this value are removed and not considered for the signature.
#' @param noisefilter Sets an upper bound for `noisefraction` cut-off below
#'   which gene expression is set to 0. Essentially gene expression above this
#'   level must be retained in the signature. Setting this higher can allow more
#'   suppression via `noisefraction` and can favour more highly expressed genes.
#' @param noisefraction Numeric value. Maximum mean log2 gene expression across
#'   cell types is calculated and values in celltypes below this fraction are
#'   set to 0. Set in conjunction with `noisefilter.` Note: if this is set too
#'   high (too close to 1), it can have a deleterious effect on deconvolution.
#' @param min_cells Numeric value specifying minimum number of cells in a
#'   subclass category. Subclass categories with fewer cells will be ignored.
#' @param remove_subclass Character vector of `subclass` levels to be removed
#'   from the analysis.
#' @param dual_mean Logical whether to calculate arithmetic mean of counts as
#'   well as mean(log2(counts +1)). This is mainly useful for simulation.
#' @param meanFUN Function for applying mean which is passed to [scmean()]. Options
#'   include `logmean` (the default) or `trimmean` which is a trimmed after
#'   excluding the top/bottom 5% of values.
#' @param postFUN Optional function applied to `genemeans` matrices after mean
#'   has been calculated. If `meanFUN` is set to `trimmean`, then `postFUN`
#'   needs to be set to `log2s`. See [scmean()].
#' @param verbose Logical whether to show messages.
#' @param sliceMem Max amount of memory in GB to allow for each subsetted count
#'   matrix object. When `scdata` is subsetted by each cell subclass, if the
#'   amount of memory would be above `sliceMem` then slicing is activated and
#'   the subsetted count matrix is divided into chunks and processed separately.
#'   This is indicated by addition of '...' in the printed timings. The limit is
#'   just under 17.2 GB (2^34 / 1e9). Above this the subsetted matrix breaches
#'   the long vector limit (>2^31 elements).
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation has increased memory requirements. See [scmean()].
#' @param ... Additional arguments passed to [scmean()].
#' @returns 
#' A list object with S3 class 'cellMarkers' containing:
#'   \item{call}{the matched call}
#'   \item{best_angle}{named list containing a matrix for each cell type with 
#'   genes in rows. Rows are ranked by lowest specificity angle for that cell 
#'   type and highest maximum expression. Columns are:
#'   `angle` the specificity angle in radians,
#'   `angle.deg` the same angle in degrees,
#'   `max` the maximum mean expression across all cell types,
#'   `rank` the rank of the mean gene expression for that cell type compared to 
#'   the other cell types}
#'   \item{group_angle}{named list of matrices similar to `best_angle`, for each 
#'   cell subclass}
#'   \item{geneset}{character vector of selected gene markers for cell types}
#'   \item{group_geneset}{character vector of selected gene markers for cell 
#'   subclasses}
#'   \item{genemeans}{matrix of mean log2+1 gene expression with genes in rows
#'   and cell types in columns}
#'   \item{genemeans_filtered}{matrix of gene expression for cell types 
#'   following noise reduction}
#'   \item{groupmeans}{matrix of mean log2+1 gene expression with genes in rows
#'   and cell subclasses in columns}
#'   \item{groupmeans_filtered}{matrix of gene expression for cell subclasses 
#'   following noise reduction}
#'   \item{cell_table}{factor encoded vector containing the groupings of the 
#'   cell types within cell subclasses, determined by which subclass contains 
#'   the maximum number of cells for each cell type}
#'   \item{spillover}{matrix of spillover values between cell types}
#'   \item{subclass_table}{contingency table of the number of cells in each
#'   subclass}
#'   \item{qqmap}{'qqmap' class object generated by [quantile_map()] containing 
#'   quantile mapping from bulk to single-cell}
#'   \item{opt}{list storing options, namely arguments `nsubclass`, `ngroup`,
#'   `expfilter`, `noisefilter`, `noisefraction`}
#'   \item{genemeans_ar}{if `dual_mean` is `TRUE`, optional matrix of arithmetic 
#'   mean, i.e. log2(mean(counts)+1)}
#'   \item{genemeans_filtered_ar}{optional matrix of arithmetic mean
#'   following noise reduction}
#' The 'cellMarkers' object is designed to be passed to [deconvolute()] to
#' deconvolute bulk RNA-Seq data. It can be updated rapidly with different
#' settings using [updateMarkers()]. Ensembl gene ids can be substituted for
#' recognisable gene symbols by applying [gene2symbol()].
#' @seealso [deconvolute()] [updateMarkers()] [gene2symbol()]
#' @author Myles Lewis
#' @importFrom matrixStats rowMaxs
#' @export
#' 
cellMarkers <- function(scdata,
                        bulkdata = NULL,
                        subclass,
                        cellgroup = NULL,
                        nsubclass = 5,
                        ngroup = 5,
                        expfilter = 0.5,
                        noisefilter = 2,
                        noisefraction = 0.25,
                        min_cells = 10,
                        remove_subclass = NULL,
                        dual_mean = FALSE,
                        meanFUN = logmean,
                        postFUN = NULL,
                        verbose = TRUE,
                        sliceMem = 16,
                        cores = 1L, ...) {
  .call <- match.call()
  if (!inherits(scdata, c("dgCMatrix", "matrix", "Seurat", "DelayedMatrix"))) {
    scdata <- as.matrix(scdata)
  }
  if (is.null(rownames(scdata))) stop("scdata is missing rownames/ gene ids")
  dimx <- dim(scdata)
  if (!is.factor(subclass)) subclass <- factor(subclass)
  if (min_cells > 0) {
    tab <- table(subclass)
    if (any(tab) < min_cells) {
      rem <- names(which(tab < min_cells))
      subclass[subclass %in% rem] <- NA
      subclass <- factor(subclass)
    }
  }
  if (!is.null(remove_subclass)) {
    subclass[subclass %in% remove_subclass] <- NA
    subclass <- factor(subclass)
  }
  if (verbose) mem_estimate(dimx, subclass, cellgroup, sliceMem, cores)
  
  ok <- TRUE
  if (!is.null(bulkdata)) {
    if (inherits(bulkdata, "data.frame")) bulkdata <- as.matrix(bulkdata)
    ok <- rownames(scdata) %in% rownames(bulkdata)
    if (verbose) message("Removing ", sum(!ok), " genes not found in bulkdata")
  }
  nsub <- nlevels(subclass)
  if (verbose) message(dimx[1], " genes, ", dimx[2], " cells, ",
                       nsub, " cell subclasses")
  if (verbose) message("Subclass analysis")
  
  nsubclass2 <- rep_len(nsubclass, nsub)
  
  if (dual_mean) {
    gm <- scmean2(scdata, subclass, meanFUN, postFUN, verbose, sliceMem, cores)
    genemeans <- gm[[1]]
    genemeans_ar <- gm[[2]]
  } else {
    genemeans <- scmean(scdata, subclass, meanFUN, postFUN, verbose, sliceMem,
                        cores, ...)
  }
  
  if (any(!ok)) {
    genemeans <- genemeans[ok, ]
    if (dual_mean) genemeans_ar <- genemeans_ar[ok, ]
    dimx[1] <- nrow(genemeans)
  }
  highexp <- rowMaxs(genemeans) > expfilter
  genemeans_filtered <- reduceNoise(genemeans[highexp, ], noisefilter,
                                    noisefraction)
  if (dual_mean) {
    genemeans_filtered_ar <- reduceNoise(genemeans_ar[highexp, ], noisefilter,
                                         noisefraction)
  }
  best_angle <- gene_angle(genemeans_filtered)
  geneset <- lapply(seq_along(best_angle), function(i) {
    rownames(best_angle[[i]])[seq_len(nsubclass2[i])]
  })
  geneset <- unique(unlist(geneset))
  
  if (!is.null(cellgroup)) {
    if (!is.factor(cellgroup)) cellgroup <- factor(cellgroup)
    # cellgroup[is.na(subclass)] <- NA
    if (min_cells > 0) {
      tab <- table(cellgroup)
      if (any(tab) < min_cells) {
        rem <- names(which(tab < min_cells))
        cellgroup[cellgroup %in% rem] <- NA
        cellgroup <- factor(cellgroup)
      }
    }
    if (verbose) message("Basic cell group analysis")
    
    # test nesting
    tab <- table(subclass, cellgroup)
    groupmeans <- scmean(scdata, cellgroup, meanFUN, postFUN, verbose,
                         sliceMem, cores, ...)
    if (any(!ok)) {
      groupmeans <- groupmeans[ok, ]
    }
    highexp <- rowMaxs(groupmeans) > expfilter
    groupmeans_filtered <- reduceNoise(groupmeans[highexp, ], noisefilter,
                                       noisefraction)
    # add extra rows
    rn <- rownames(groupmeans_filtered)
    miss <- rn[!rn %in% rownames(genemeans_filtered)]
    if (length(miss) > 0) {
      extra <- reduceNoise(genemeans[miss, , drop = FALSE], noisefilter,
                           noisefraction)
      genemeans_filtered <- rbind(genemeans_filtered, extra)
    }
    group_angle <- gene_angle(groupmeans_filtered)
    ngroup2 <- rep_len(ngroup, length(group_angle))
    group_geneset <- lapply(seq_along(group_angle), function(i) {
      rownames(group_angle[[i]])[seq_len(ngroup2[i])]
    })
    group_geneset <- unique(unlist(group_geneset))
    cell_table <- apply(tab, 1, function(i) names(which.max(i)))
    cell_table <- factor(cell_table, levels = unique(cell_table))
  } else {
    group_geneset <- group_angle <- groupmeans <- groupmeans_filtered <- NULL
    cell_table <- NULL
  }
  
  # determine spillover
  gene_sig <- genemeans_filtered[geneset, ]
  m_itself <- dotprod(gene_sig, gene_sig, equal_weight = FALSE)
  
  # map bulk to sc
  qqmap <- if (!is.null(bulkdata)) {
    message("Quantile map bulk to sc, ", appendLF = FALSE)
    quantile_map(log2(bulkdata +1), genemeans, remove_zeros = TRUE)
  } else NULL
  
  out <- list(call = .call,
              best_angle = best_angle,
              group_angle = group_angle,
              geneset = geneset,
              group_geneset = group_geneset,
              genemeans = genemeans,
              genemeans_filtered = genemeans_filtered,
              groupmeans = groupmeans,
              groupmeans_filtered = groupmeans_filtered,
              cell_table = cell_table,
              spillover = m_itself,
              subclass_table = table(subclass, dnn = NULL),
              qqmap = qqmap,
              opt = list(nsubclass = nsubclass,
                         ngroup = ngroup,
                         expfilter = expfilter,
                         noisefilter = noisefilter,
                         noisefraction = noisefraction))
  if (dual_mean) {
    out$genemeans_ar <- genemeans_ar
    out$genemeans_filtered_ar <- genemeans_filtered_ar
  }
  class(out) <- "cellMarkers"
  out
}


#' @export
summary.cellMarkers <- function(object, ...) {
  cat("scRNA-Seq dataset:", nrow(object$genemeans), "genes\n")
  cat("Cell subclasses:", ncol(object$genemeans), "\n")
  cat("Cell subclass signature:", length(object$geneset), "genes\n")
  cat("Cell groups:", ncol(object$groupmeans), "\n")
  cat("Cell group signature:", length(object$group_geneset), "genes\n")
  mmeth <- if (is.null(object$call[["meanFUN"]])) {
    "logmean"
  } else as.character(object$call[["meanFUN"]])
  cat("Mean method:", mmeth, "\n")
  for (i in 1:5) {
    cat(paste0(names(object$opt)[i], ": ", object$opt[i], "\n"))
  }
  cat("Subclass cell totals:\n")
  print(object$subclass_table)
}


# estimate memory requirement
mem_estimate <- function(dimx, subclass, cellgroup, sliceMem, cores) {
  tab <- table(subclass)
  tab2 <- table(cellgroup)
  dimmax <- as.numeric(max(c(tab, tab2), na.rm = TRUE)) * dimx[1]
  mem <- structure(dimmax * 8, class = "object_size")
  if (mem > 1e9)
    message("Max subclass/group memory ", format(mem, units = "GB"))
  if (mem > sliceMem * 1e9) {
    message("Slicing above ", sliceMem, " Gb")
  } else message("No slicing")
  
  nsubcl <- sort(tab, decreasing = TRUE)[1:cores]
  ngroup <- sort(tab2, decreasing = TRUE)[1:cores]
  msubcl <- nsubcl * as.numeric(dimx[1]) / 1.25e8
  mgroup <- ngroup * as.numeric(dimx[1]) / 1.25e8
  msubcl[msubcl > sliceMem] <- sliceMem
  mgroup[mgroup > sliceMem] <- sliceMem
  mem <- max(c(sum(msubcl), sum(mgroup)))
  message("Estimated memory requirement ", format(mem, digits = 2), " Gb")
}
