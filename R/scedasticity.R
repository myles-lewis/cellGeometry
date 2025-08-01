
#' Residuals plot
#' 
#' Plots residuals from a deconvolution result object against bulk gene
#' expression (on semi-log axis). Normal residuals, weighted residuals or
#' Studentized residuals can be visualised to check for heteroscedasticity and
#' genes with extreme errors.
#' 
#' @param fit 'deconv' class deconvolution object
#' @param test bulk gene expression matrix assumed to be in raw counts
#' @param type Specifies type of residuals to be plotted
#' @param ... Optional arguments passed to [plot()]
#' @returns Produces a scatter plot in base graphics. Returns invisibly a
#'   dataframe of the coordinates of the points.
#' @export
plot_residuals <- function(fit, test, type = c("reg", "student", "weight"),
                           ...) {
  if (!inherits(fit, "deconv")) stop("not a 'deconv' class object")
  nm <- fit$call$test
  .call <- match.call()
  if (nm != .call$test) {
    message("`", .call$test, "` does not match the bulk dataset `", nm,
            "` used in `", .call$fit, "`")
  }
  type <- match.arg(type)
  geneset <- rownames(fit$subclass$residuals)
  ge <- test[geneset, ]
  res <- fit$subclass$residuals
  res <- switch(type, student = rstudent(fit),
                weight = res * fit$subclass$weights,
                res)
  ylab <- switch(type, student = "Studentized residuals",
                 weight = "Weighted residuals",
                 "Residuals")
  new.args <- list(...)
  args <- list(x = ge, y = res, log = "x",
               xlim = c(1, max(ge, na.rm = TRUE)),
               cex = 0.7, col = adjustcolor("black", 0.2), pch = 16, bty = "l",
               xlab = "Bulk gene expression", ylab = ylab)
  if (length(new.args)) args[names(new.args)] <- new.args     
  do.call("plot", args) |>
    suppressWarnings()
  abline(0, 0, col = "red")
  invisible(data.frame(expr = as.vector(ge), res = as.vector(res),
                       gene = rownames(res)))
}
