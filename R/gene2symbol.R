
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
#' @returns If `x` is a vector, a vector of symbols is returned. If no symbol is
#'   no available the ensembl id is left untouched. If `x` is a 'cellMarkers'
#'   class object, a 'cellMarkers' object is returned with rownames in the
#'   results elements converted to gene symbols.
#' @export

gene2symbol <- function(x, ensdb) {
  if (inherits(x, "cellMarkers")) {
    rownames(x$genemeans) <- convertsymbol(rownames(x$genemeans), ensdb)
    rownames(x$genemeans_filtered) <- convertsymbol(rownames(x$genemeans_filtered), ensdb)
    rownames(x$groupmeans) <- convertsymbol(rownames(x$groupmeans), ensdb)
    rownames(x$groupmeans_filtered) <- convertsymbol(rownames(x$groupmeans_filtered), ensdb)
    x$geneset <- convertsymbol(x$geneset, ensdb)
    x$group_geneset<- convertsymbol(x$group_geneset, ensdb)
    return(x)
  }
  convertsymbol(x, ensdb)
}


#' @importFrom ensembldb select
#' 
convertsymbol <- function(x, ensdb) {
  geneid <- ensembldb::select(ensdb, keys = x, keytype = "GENEID",
                              columns = c("GENEID", "SYMBOL"))
  ok <- geneid$SYMBOL != "" & !is.na(geneid$SYMBOL)
  x[ok] <- geneid$SYMBOL[ok]
  attr(x, "ok") <- ok
  x
}
