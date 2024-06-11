
#' Update cellMarkers object
#'
#' Updates a cellMarkers gene signature object without having to rerun
#' calculation of gene means which can be slow.
#'
#' @param object A cellMarkers class object. Either `object` or `genemeans` must
#'   be specified.
#' @param genemeans A matrix of mean gene expression with genes in rows and cell
#'   subclasses in columns.
#' @param groupmeans Optional matrix of mean gene expression for overarching
#'   main cell groups (genes in rows, cell groups in columns).
#' @param add_genes Character vector of gene markers to add manually to the cell
#'   subclass gene signature.
#' @param add_groupgenes Character vector of gene markers to add manually to the
#'   cell group gene signature.
#' @param nsubclass Number of genes to select for each single cell subclass.
#' @param ngroup Number of genes to select for each cell group.
#' @param expfilter Genes whose maximum mean expression on log2 scale per cell
#'   type are below this value are removed.
#' @param noisefilter Numeric value below which mean log2 gene expression is
#'   reduced to 0.
#' @param noisefraction Numeric value. Maximum mean log2 gene expression across
#'   cell types is calculated and values in celltypes below this fraction are
#'   set to 0.
#' @param verbose Logical whether to show messages.
#' @returns Returns list object containing `best_angle`, a list of genes ranked
#' by lowest angle and highest maximum expression in a cell type; `genemeans`,
#' matrix of mean log2+1 gene expression with genes in rows and cell types in
#' columns; `genemeans_filters`, matrix of gene expression following noise
#' reduction.
#' @export

updateMarkers <- function(object = NULL,
                          genemeans = NULL,
                          groupmeans = NULL,
                          add_genes = NULL,
                          add_groupgenes = NULL,
                          nsubclass = 5,
                          ngroup = 5,
                          expfilter = 1,
                          noisefilter = 1.5,
                          noisefraction = 0.25,
                          verbose = TRUE) {
  .call <- match.call()
  
  if (is.null(object) && is.null(genemeans))
    stop("Either a cellMarkers object or genemeans must be supplied")
  
  if (is.null(genemeans)) genemeans <- object$genemeans
  if (any(duplicated(rownames(genemeans))))
    stop("Duplicated rownames in genemeans")
  
  if (verbose) message("Subclass analysis")
  highexp <- rowMaxs(genemeans) > expfilter |
    rownames(genemeans) %in% c(add_genes, add_groupgenes)
  genemeans_filtered <- reduceNoise(genemeans[highexp, ], noisefilter,
                                    noisefraction)
  best_angle <- gene_angle(genemeans_filtered)
  geneset <- lapply(best_angle, function(i) rownames(i)[1:nsubclass])
  geneset <- unique(c(unlist(geneset), add_genes))
  
  if (is.null(groupmeans)) groupmeans <- object$groupmeans
  
  if (!is.null(groupmeans)) {
    if (any(duplicated(rownames(groupmeans))))
      stop("Duplicated rownames in groupmeans")
    if (verbose) message("Basic cell group analysis")
    
    highexp <- rowMaxs(groupmeans) > expfilter
    groupmeans_filtered <- reduceNoise(groupmeans[highexp, ], noisefilter,
                                       noisefraction)
    group_angle <- gene_angle(groupmeans_filtered)
    group_geneset <- lapply(group_angle, function(i) rownames(i)[1:ngroup])
    group_geneset <- unique(c(unlist(group_geneset), add_groupgenes))
    cell_table <- object$cell_table
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
