
#' Update cellMarkers object
#'
#' Updates a 'cellMarkers' gene signature object with new settings without
#' having to rerun calculation of gene means, which can be slow.
#'
#' @param object A 'cellMarkers' class object. Either `object` or `genemeans`
#'   must be specified.
#' @param genemeans A matrix of mean gene expression with genes in rows and cell
#'   subclasses in columns.
#' @param groupmeans Optional matrix of mean gene expression for overarching
#'   main cell groups (genes in rows, cell groups in columns).
#' @param add_gene Character vector of gene markers to add manually to the cell
#'   subclass gene signature.
#' @param add_groupgene Character vector of gene markers to add manually to the
#'   cell group gene signature.
#' @param remove_gene Character vector of gene markers to manually remove from
#'   the cell subclass gene signature.
#' @param remove_groupgene Character vector of gene markers to manually remove
#'   to the cell group gene signature.
#' @param remove_subclass Character vector of cell subclasses to remove.
#' @param remove_group Optional character vector of cell groups to remove.
#' @param bulkdata Optional data matrix containing bulk RNA-Seq data with genes
#'   in rows. This matrix is only used for its rownames, to ensure that cell
#'   markers are selected from genes in the bulk dataset.
#' @param nsubclass Number of genes to select for each single cell subclass.
#'   Either a single number or a vector with the number of genes for each
#'   subclass.
#' @param ngroup Number of genes to select for each cell group.
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
#' @param verbose Logical whether to show messages.
#' @returns A list object of S3 class 'cellMarkers'. See [cellMarkers()] for
#'   details. If [gene2symbol()] has been called, an extra list element `symbol`
#'   will be present. The list element `update` stores the call to 
#'   `updateMarkers()`.
#' @seealso [cellMarkers()] [gene2symbol()]
#' @author Myles Lewis
#' @export

updateMarkers <- function(object = NULL,
                          genemeans = NULL,
                          groupmeans = NULL,
                          add_gene = NULL,
                          add_groupgene = NULL,
                          remove_gene = NULL,
                          remove_groupgene = NULL,
                          remove_subclass = NULL,
                          remove_group = NULL,
                          bulkdata = NULL,
                          nsubclass = object$opt$nsubclass,
                          ngroup = object$opt$ngroup,
                          expfilter = object$opt$expfilter,
                          noisefilter = object$opt$noisefilter,
                          noisefraction = object$opt$noisefraction,
                          verbose = TRUE) {
  .call <- match.call()
  
  if (is.null(object) && is.null(genemeans))
    stop("Either a cellMarkers object or genemeans must be supplied")
  
  if (is.null(genemeans)) genemeans <- object$genemeans
  if (is.null(groupmeans)) groupmeans <- object$groupmeans
  if (any(duplicated(rownames(genemeans))))
    stop("Duplicated rownames in genemeans")
  
  ok <- TRUE
  if (!is.null(bulkdata)) {
    ok <- rownames(genemeans) %in% rownames(bulkdata)
    if (any(!ok)) {
      if (verbose) message(sum(ok), " genes overlap with bulkdata")
    }
  }
  
  # remove subclass or group
  if (!is.null(remove_subclass) | !is.null(remove_group)) {
    if (any(!remove_subclass %in% colnames(genemeans))) stop("cannot remove subclass")
    if (any(!remove_group %in% colnames(groupmeans))) stop("cannot remove group")
    
    subcl <- !colnames(genemeans) %in% remove_subclass &
      !object$cell_table %in% remove_group
    genemeans <- genemeans[, subcl]
    reduced_tab <- object$cell_table[subcl]
    old_levels <- levels(object$cell_table)
    object$cell_table <- droplevels(reduced_tab)
    object$subclass_table <- object$subclass_table[subcl]
    gone <- setdiff(old_levels, levels(object$cell_table))
    remove_group <- c(remove_group, gone)
    grp <- !colnames(object$groupmeans) %in% remove_group
    groupmeans <- if (sum(grp) > 1) groupmeans[, grp] else NULL
    object$genemeans_ar <- object$genemeans_ar[, subcl]
  }
  
  # subclass analysis
  if (verbose) message("Subclass analysis")
  nsub <- length(object$subclass_table)
  nsubclass2 <- rep_len(nsubclass, nsub)
  highexp <- ok & rowMaxs(genemeans) > expfilter |
    rownames(genemeans) %in% add_gene
  genemeans_filtered <- reduceNoise(genemeans[highexp, ], noisefilter,
                                    noisefraction)
  best_angle <- gene_angle(genemeans_filtered)
  geneset <- lapply(seq_along(best_angle), function(i) {
    rownames(best_angle[[i]])[seq_len(nsubclass2[i])]
  })
  geneset <- unique(c(unlist(geneset), add_gene))
  if (!is.null(remove_gene)) geneset <- geneset[!geneset %in% remove_gene]
  geneset <- geneset[!is.na(geneset)]
  if (verbose) message(length(geneset), " marker genes")
  
  # group analysis
  if (!is.null(groupmeans)) {
    if (any(duplicated(rownames(groupmeans))))
      stop("Duplicated rownames in groupmeans")
    if (verbose) message("Basic cell group analysis")
    
    highexp <- ok & rowMaxs(groupmeans) > expfilter |
      rownames(groupmeans) %in% add_groupgene
    groupmeans_filtered <- reduceNoise(groupmeans[highexp, ], noisefilter,
                                       noisefraction)
    # add extra rows
    rn <- rownames(groupmeans_filtered)
    miss <- rn[!rn %in% rownames(genemeans_filtered)]
    if (length(miss) > 0) {
      extra <- reduceNoise(genemeans[miss, , drop = FALSE], noisefilter, noisefraction)
      genemeans_filtered <- rbind(genemeans_filtered, extra)
    }
    group_angle <- gene_angle(groupmeans_filtered)
    ngroup2 <- rep_len(ngroup, length(group_angle))
    group_geneset <- lapply(seq_along(group_angle), function(i) {
      rownames(group_angle[[i]])[seq_len(ngroup2[i])]
    })
    group_geneset <- unique(c(unlist(group_geneset), add_groupgene))
    if (!is.null(remove_groupgene)) {
      group_geneset <- group_geneset[!group_geneset %in% remove_groupgene]
    }
    cell_table <- object$cell_table
  } else {
    group_geneset <- group_angle <- groupmeans <- groupmeans_filtered <- NULL
    cell_table <- NULL
  }
  
  # determine spillover
  gene_sig <- genemeans_filtered[geneset, ]
  m_itself <- dotprod(gene_sig, gene_sig)
  
  out <- list(call = object$call,
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
              subclass_table = object$subclass_table,
              opt = list(nsubclass = nsubclass,
                         ngroup = ngroup,
                         expfilter = expfilter,
                         noisefilter = noisefilter,
                         noisefraction = noisefraction),
              update = .call)
  if (!is.null(object$symbol)) out$symbol <- object$symbol
  if (!is.null(object$qqmerge)) out$qqmerge <- object$qqmerge
  
  if (!is.null(object$genemeans_ar)) {
    # dual mean
    genemeans_ar <- object$genemeans_ar
    highexp <- ok & rowMaxs(genemeans_ar) > expfilter |
      rownames(genemeans_ar) %in% add_gene
    out$genemeans_filtered_ar <- reduceNoise(genemeans_ar[highexp, ],
                                             noisefilter, noisefraction)
    out$genemeans_ar <- genemeans_ar
  }
  class(out) <- "cellMarkers"
  out
}
