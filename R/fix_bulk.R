
#' Fix in missing genes in bulk RNA-Seq matrix 
#'
#' Fills in missing genes in a bulk RNA-Seq matrix based on the gene signature
#' of a 'cellMarkers' objects. Signature is taken from both the subclass gene
#' set and group gene set.
#' 
#' @param mk object of class 'cellMarkers'. See [cellMarkers()].
#' @param bulk matrix of bulk RNA-Seq
#' @returns Expanded bulk matrix with extra rows for missing genes, filled with
#'   zeros.
#' @export

fix_bulk <- function(mk, bulk) {
  if (!inherits(mk, "cellMarkers")) stop("Not a 'cellMarkers' class object")
  genes <- unique(c(mk$geneset, mk$group_geneset))
  ok <- genes %in% rownames(bulk)
  if (all(ok)) return(bulk)
  pc <- format(sum(!ok) / length(genes) *100, digits = 2)
  message(sum(!ok), "/", length(genes), " (", pc, 
          "%) signature genes missing from bulk")
  miss <- genes[!ok]
  extra <- matrix(0, nrow = length(miss), ncol = ncol(bulk))
  rownames(extra) <- miss
  rbind(as.matrix(bulk), extra)
}
