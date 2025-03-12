


subMarkers <- function(object, subset, nsubclass = NULL,
                       expfilter = NULL, noisefilter = NULL,
                       noisefraction = NULL) {
  if (!inherits(object, "cellMarkers")) stop("not a 'cellMarkers' object")
  if (is.null(nsubclass)) nsubclass <- object$opt$nsubclass
  if (is.null(expfilter)) expfilter <- object$opt$expfilter
  if (is.null(noisefilter)) noisefilter <- object$opt$noisefilter
  if (is.null(noisefraction)) noisefraction <- object$opt$noisefraction
  
  if (!subset %in% levels(object$cell_table)) stop("subset not found")
  w <- object$cell_table %in% subset
  
  genemeans <- object$genemeans_filtered[, w]
  highexp <- rowMaxs(genemeans) > expfilter
  genemeans_filtered <- reduceNoise(genemeans[highexp, ], noisefilter,
                                    noisefraction)
  best_angle <- gene_angle(genemeans_filtered)
  geneset <- lapply(seq_along(best_angle), function(i) {
    rownames(best_angle[[i]])[seq_len(nsubclass)]
  })
  geneset <- unique(unlist(geneset))
  list(geneset = geneset, best_angle = best_angle,
       genemeans_filtered = genemeans_filtered)
}
