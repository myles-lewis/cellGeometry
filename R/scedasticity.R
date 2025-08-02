
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
#' @param show_outliers Logical whether to show any remaining outlying extreme
#'   genes in red
#' @param ... Optional arguments passed to [plot()]
#' @returns Produces a scatter plot in base graphics. Returns invisibly a
#'   dataframe of the coordinates of the points.
#' @export
plot_residuals <- function(fit, test, type = c("reg", "student", "weight"),
                           show_outliers = TRUE,
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
  col <- adjustcolor("black", 0.2)
  dat <- data.frame(expr = as.vector(ge), res = as.vector(res),
                    gene = rownames(res))
  if (show_outliers) {
    outlier_method <- fit$call$outlier_method %||% "var.e"
    outlier_cutoff <- fit$call$outlier_cutoff %||% switch(outlier_method, var.e = 4,
                                                          cooks = 1, rstudent = 10)
    outlier_quantile <- fit$call$outlier_quantile %||% 0.9
    count_space <- fit$call$count_space %||% TRUE
    metric <- outlier_metric(fit$subclass, outlier_method, outlier_quantile, count_space)
    outlier <- metric > outlier_cutoff
    nm <- names(metric[outlier])
    dat$outlier <- FALSE
    dat$outlier[dat$gene %in% nm] <- TRUE
    dat <- dat[order(dat$outlier), ]
    col <- rep(col, nrow(dat))
    col[dat$outlier] <- adjustcolor("red", 0.7)
  }
  new.args <- list(...)
  args <- list(x = dat$expr, y = dat$res, log = "x",
               xlim = c(pmax(1, min(ge, na.rm = TRUE)), max(ge, na.rm = TRUE)),
               cex = 0.7, col = col, pch = 16, bty = "l",
               xlab = "Bulk gene expression", ylab = ylab)
  if (length(new.args)) args[names(new.args)] <- new.args     
  do.call("plot", args) |>
    suppressWarnings()
  abline(0, 0, col = "blue")
  invisible(dat)
}
