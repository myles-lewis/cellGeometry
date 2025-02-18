
#' Converts ensembl gene ids to symbols 
#'
#' Uses a loaded ensembl database to convert ensembl gene ids to symbol. If a
#' vector is provided, a vector of symbols is returned. If a cellMarkers object
#' is provided, the rownames in the genemeans, genemeans_filtered, groupmeans
#' and groupmeans_filtered elements are changed to symbol and the cellMarkers
#' object is returned.
#' 
#' @param x Either a vector of ensembl gene ids to convert or a 'cellMarkers'
#'   class object.
#' @param ensdb An ensembl database object loaded via the `AnnotationHub`
#'   bioconductor package.
#' @param dups Character vector specifying action for duplicated gene symbols.
#'   "omit" means that duplicated gene symbols are not replaced, but left as
#'   ensembl gene ids. "pass" means that all gene ids are replaced where
#'   possible even if that leads to duplicates. Duplicates can cause problems
#'   with rownames and [updateMarkers()] in particular.
#' @returns If `x` is a vector, a vector of symbols is returned. If no symbol is
#'   available for particular ensembl id, the id is left untouched. If `x` is a
#'   'cellMarkers' class object, a 'cellMarkers' object is returned with
#'   rownames in the results elements and genesets converted to gene symbols,
#'   and an extra element `symbol` containing a named vector of converted genes.
#' @seealso [cellMarkers()]
#' @export

gene2symbol <- function(x, ensdb, dups = c("omit", "pass")) {
  dups <- match.arg(dups)
  if (inherits(x, "cellMarkers")) {
    old_rn <- rownames(x$genemeans)
    rn <- rownames(x$genemeans) <- convertsymbol(old_rn, ensdb, dups)
    names(rn) <- old_rn
    rownames(x$groupmeans) <- rn
    rownames(x$genemeans_filtered) <- rn[rownames(x$genemeans_filtered)]
    rownames(x$groupmeans_filtered) <- rn[rownames(x$groupmeans_filtered)]
    for (i in seq_along(x$best_angle)) {
      rownames(x$best_angle[[i]]) <- rn[rownames(x$best_angle[[i]])]
    }
    for (i in seq_along(x$group_angle)) {
      rownames(x$group_angle[[i]]) <- rn[rownames(x$group_angle[[i]])]
    }
    if (!is.null(x$genemeans_ar)) {
      rownames(x$genemeans_ar) <- rn
      rownames(x$genemeans_filtered_ar) <- rn[rownames(x$genemeans_filtered_ar)]
    }
    
    x$geneset <- rn[x$geneset]
    x$group_geneset<- rn[x$group_geneset]
    x$symbol <- rn
    return(x)
  }
  convertsymbol(x, ensdb)
}


#' @importFrom ensembldb select
#' 
convertsymbol <- function(x, ensdb, dups = "omit") {
  geneid <- ensembldb::select(ensdb, keys = x, keytype = "GENEID",
                              columns = c("GENEID", "SYMBOL"))
  ok <- which(geneid$SYMBOL != "" & !is.na(geneid$SYMBOL))
  newx <- geneid$SYMBOL[ok]
  if (dups == "omit") {
    dups <- duplicated(newx)
    x[ok[!dups]] <- newx[!dups]
  } else {
    x[ok] <- newx
  }
  x
}
