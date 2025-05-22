

fix_bulk <- function(mk, bulk) {
  if (!inherits(mk, "cellMarkers")) stop("Not a 'cellMarkers' class object")
  genes <- unique(c(mk$geneset, mk$group_geneset))
  ok <- genes %in% rownames(bulk)
  if (all(ok)) return(bulk)
  pc <- format(sum(!ok) / length(genes), digits = 2)
  message(sum(!ok), "/", length(genes), " (", pc, 
          "%) signature genes missing from bulk")
  miss <- genes[!ok]
  extra <- matrix(0, nrow = length(miss), ncol = ncol(bulk))
  rownames(extra) <- miss
  rbind(bulk, extra)
}
