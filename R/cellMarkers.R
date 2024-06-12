
#' Identify cell markers
#' 
#' Uses geometric method based on vector dot product to identify genes which are
#' the best markers for individual cell types.
#' 
#' @param scdata Single-cell data matrix with genes in rows and cells in
#'   columns.
#' @param bulkdata Optional data matrix containing bulk RNA-Seq data with genes
#'   in rows. This matrix is only used for its rownames, to ensure that cell
#'   markers are selected from genes in the bulk dataset.
#' @param subclass Vector of cell subclasses matching the columns in `scdata`
#' @param cellgroup Optional grouping vector of major cell types matching the
#'   columns in `scdata`. `subclass` is assumed to contain subclasses which are
#'   subsets within `cellgroup` overarching classes.
#' @param nsubclass Number of genes to select for each single cell subclass.
#' @param ngroup Number of genes to select for each cell group. 
#' @param expfilter Genes whose maximum mean expression on log2 scale per cell
#'   type are below this value are removed.
#' @param noisefilter Numeric value below which mean log2 gene expression is
#'   reduced to 0.
#' @param noisefraction Numeric value. Maximum mean log2 gene expression across
#'   cell types is calculated and values in celltypes below this fraction are
#'   set to 0.
#' @param big Logical whether to invoke matrix slicing to handle big matrices.
#' @param verbose Logical whether to show messages.
#' @returns 
#' Returns list object containing `best_angle`, a list of genes ranked by lowest
#' angle and highest maximum expression in a cell type; `genemeans`, matrix of
#' mean log2+1 gene expression with genes in rows and cell types in columns;
#' `genemeans_filters`, matrix of gene expression following noise reduction.
#' @importFrom matrixStats rowMaxs
#' @export
#' 
cellMarkers <- function(scdata,
                        bulkdata = NULL,
                        subclass,
                        cellgroup = NULL,
                        nsubclass = 5,
                        ngroup = 5,
                        expfilter = 1,
                        noisefilter = 1.5,
                        noisefraction = 0.25,
                        big = NULL,
                        verbose = TRUE) {
  .call <- match.call()
  if (!inherits(scdata, c("dgCMatrix", "matrix", "Seurat"))) scdata <- as.matrix(scdata)
  dimx <- dim(scdata)
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  if (!is.factor(subclass)) subclass <- factor(subclass)
  
  if (!is.null(bulkdata)) {
    ok <- rownames(scdata) %in% rownames(bulkdata)
    if (any(!ok) & (is.null(big) || !big)) {
      if (verbose) message("Removing ", sum(!ok), " genes not found in bulkdata")
      scdata <- scdata[ok, ]
      dimx <- dim(scdata)
    }
  }
  u <- unique(subclass)
  nsub <- length(u[!is.na(u)])
  if (verbose) message(dimx[1], " genes, ", dimx[2], " cells, ",
                       nsub, " cell subclasses")
  if (verbose) message("Subclass analysis")
  genemeans <- scmean(scdata, subclass, big, verbose)
  highexp <- rowMaxs(genemeans) > expfilter
  genemeans_filtered <- reduceNoise(genemeans[highexp, ], noisefilter,
                                    noisefraction)
  best_angle <- gene_angle(genemeans_filtered)
  geneset <- lapply(best_angle, function(i) rownames(i)[1:nsubclass])
  geneset <- unique(unlist(geneset))
  
  if (!is.null(cellgroup)) {
    if (verbose) message("Basic cell group analysis")
    if (!is.factor(cellgroup)) cellgroup <- factor(cellgroup)
    # test nesting
    tab <- table(subclass, cellgroup)
    groupmeans <- scmean(scdata, cellgroup, big, verbose)
    highexp <- rowMaxs(groupmeans) > expfilter
    groupmeans_filtered <- reduceNoise(groupmeans[highexp, ], noisefilter,
                                       noisefraction)
    group_angle <- gene_angle(groupmeans_filtered)
    group_geneset <- lapply(group_angle, function(i) rownames(i)[1:ngroup])
    group_geneset <- unique(unlist(group_geneset))
    cell_table <- apply(tab, 1, function(i) names(which.max(i)))
    cell_table <- factor(cell_table, levels = unique(cell_table))
  } else {
    group_geneset <- group_angle <- groupmeans <- groupmeans_filtered <- NULL
    cell_table <- NULL
  }
  
  # determine spillover
  gene_sig <- genemeans_filtered[geneset, ]
  m_itself <- dotprod(gene_sig, gene_sig, equalWeight = FALSE)
  
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
              spillover = m_itself)
  class(out) <- "cellMarkers"
  out
}


scmean <- function(x, celltype, big = NULL, verbose = TRUE) {
  start0 <- Sys.time()
  if (!is.factor(celltype)) celltype <- factor(celltype)
  ok <- !is.na(celltype)
  dimx <- dim(x)
  if (dimx[2] != length(celltype)) stop("Incompatible dimensions")
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  if (is.null(big) || !big) {
    # small matrix
    genemeans <- vapply(levels(celltype), function(i) {
      logmean(x[, which(celltype==i & ok)])
    }, numeric(dimx[1]))
    return(genemeans)
  }
  # large matrix
  s <- sliceIndex(dimx[1])
  genemeans <- vapply(levels(celltype), function(i) {
    if (verbose) cat(i, " ")
    start <- Sys.time()
    c_index <- which(celltype == i & ok)
    out <- lapply(s, function(j) {
      logmean(as.matrix(x[j, c_index])) |> suppressWarnings()
    })
    if (verbose) timer(start)
    unlist(out)
  }, numeric(dimx[1]))
  
  if (verbose) timer(start0, "Duration")
  genemeans
}

logmean <- function(x) rowMeans(log2(x +1))

sliceIndex <- function(nx, sliceSize = 2000) {
  if (is.null(sliceSize)) sliceSize <- nx
  s <- ceiling(nx / sliceSize)
  excess <- nx %% sliceSize
  lapply(seq_len(s), function(i) {
    if (i==s && excess != 0) return(seq_len(excess) + sliceSize * (i-1L))
    seq_len(sliceSize) + sliceSize * (i-1L)
  })
}

timer <- function(start, msg = NULL) {
  end <- Sys.time()
  if (is.null(msg)) {
    cat(paste0("(", format(end - start, digits = 3), ")\n"))
  } else {
    cat(msg, format(end - start, digits = 3), "\n")
  }
}
