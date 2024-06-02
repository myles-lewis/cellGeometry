
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
                        ngroup = 0,
                        expfilter = 1,
                        noisefilter = 1.5,
                        noisefraction = 0.25) {
  scdata <- as.matrix(scdata)
  if (!is.factor(subclass)) subclass <- factor(subclass)
  
  if (!is.null(bulkdata)) {
    ok <- rownames(scdata) %in% rownames(bulkdata)
    if (any(!ok)) {
      message("Removing ", sum(!ok), " genes not found in bulkdata")
      scdata <- scdata[ok, ]
    }
  }
  u <- unique(subclass)
  nsub <- length(u[!is.na(u)])
  message(nrow(scdata), " genes, ", ncol(scdata), " cells, ",
          nsub, " cell subclasses")
  message("Subclass analysis")
  genemeans <- scmean(scdata, subclass)
  highexp <- rowMaxs(genemeans) > expfilter
  genemeans_filtered <- reduceNoise(genemeans[highexp, ], noisefilter,
                                    noisefraction)
  best_angle <- gene_angle(genemeans_filtered)
  geneset <- lapply(best_angle, function(i) rownames(i)[1:nsubclass])
  geneset <- unique(unlist(geneset))
  
  if (!is.null(cellgroup)) {
    message("Basic cell group analysis")
    if (!is.factor(cellgroup)) cellgroup <- factor(cellgroup)
    # test nesting
    tab <- table(subclass, cellgroup)
    tab <- tab > 0L
    if (any(rowSums(tab) != 1)) stop("subclass is not nested in cellgroup")
    
    groupmeans <- scmean(scdata, cellgroup)
    highexp <- rowMaxs(groupmeans) > expfilter
    groupmeans_filtered <- reduceNoise(groupmeans[highexp, ], noisefilter,
                                       noisefraction)
    group_angle <- gene_angle(groupmeans_filtered)
    group_geneset <- lapply(group_angle, function(i) rownames(i)[1:ngroup])
    group_geneset <- unique(unlist(group_geneset))
    cell_table <- apply(tab, 1, function(i) names(which(i)))
    cell_table <- factor(cell_table, levels = unique(cell_table))
  } else {
    group_geneset <- group_angle <- groupmeans <- groupmeans_filtered <- NULL
    cell_table <- NULL
  }
  
  # determine spillover
  gene_sig <- genemeans_filtered[geneset, ]
  m_itself <- dotprod(gene_sig, gene_sig, equalWeight = FALSE)
  
  out <- list(best_angle = best_angle,
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


scmean <- function(scdata, celltype) {
  scdata <- as.matrix(scdata)
  if (!is.factor(celltype)) celltype <- factor(celltype)
  ok <- !is.na(celltype)
  celltype <- celltype[ok]
  logcounts <- as.matrix(log2(scdata[, ok] +1))
  genemeans <- vapply(levels(celltype), function(i) rowMeans(logcounts[, celltype==i]),
                      numeric(nrow(logcounts)))
  genemeans
}


# Vector based best marker selection
gene_angle <- function(genemeans) {
  genemeans_scaled <- scaleSphere(genemeans)
  genemeans_angle <- acos(genemeans_scaled)
  genemeans_max <- rowMaxs(genemeans)
  best_angle <- lapply(colnames(genemeans_angle), function(i) {
    df <- data.frame(angle=genemeans_angle[, i], max=genemeans_max)
    df[with(df, order(angle, -max)), ]
  })
  names(best_angle) <- colnames(genemeans_angle)
  best_angle
}

# unit hypersphere scaling applied to genes
# genes in rows, celltypes in columns
scaleSphere <- function(cellmat) {
  vecLength <- sqrt(rowSums(cellmat^2))
  cellmat / vecLength
}


# hypersphere scaling of genes adjusted for total read count across cell subclasses
adjScaleGeneMatrix <- function(gene_sign, celltotals, meandepth) {
  sig_unlog <- 2^gene_sign - 1
  sig_scaleto_bulkdepth <- t(t(sig_unlog) / celltotals * meandepth)
  sig_scaled <- scaleSphere(sig_scaleto_bulkdepth) * 2^10
  log_sig <- log2(sig_scaled + 1)
  log_sig
}
