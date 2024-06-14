
#' Vector based best marker selection
#' 
#' Core function which takes a matrix of mean gene expression (assumed to be
#' log2 transformed to be more Gaussian). Mean gene expression per gene is
#' scaled to a unit hypersphere assuming each gene represents a vector in space
#' with dimensions representing each cell subclass/group.
#' 
#' @param genemeans matrix of mean gene expression with genes in rows and
#'   celltypes, tissues or subclasses in columns.
#' @returns a list whose length is the number of columns in genemeans, with each
#'   element containing a dataframe with genes in rows, sorted by best marker
#'   status as determined by minimum vector angle and highest maximum gene
#'   expression per celltype/tissue.
#' @export
#' 
gene_angle <- function(genemeans) {
  genemeans_scaled <- scaleSphere(genemeans)
  genemeans_angle <- acos(genemeans_scaled)
  genemeans_max <- rowMaxs(genemeans)
  gene_rank <- t(apply(-genemeans, 1, rank))
  best_angle <- lapply(colnames(genemeans_angle), function(i) {
    df <- data.frame(angle = genemeans_angle[, i],
                     angle.deg = genemeans_angle[, i] * 180 / pi,
                     max = genemeans_max,
                     rank = gene_rank[, i])
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
