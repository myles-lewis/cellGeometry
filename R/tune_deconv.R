
#' Tune deconvolution parameters
#' 
#' Tests a tuning grid of a deconvolution parameter for either [updateMarkers()]
#' (e.g. `expfilter` or `nsubclass`) or [deconvolute()] (e.g. `comp_amount`).
#' 
#' @param cm cellMarkers class object
#' @param test matrix of bulk RNA-Seq to be deconvoluted. Passed to [deconvolute()].
#' @param samples matrix of cell amounts with subclasses in columns and samples
#'   in rows.
#' @param grid Named list of vectors for the tuning grid similar to
#'   [expand.grid()]. Names represent the parameter to be tuned which must be an
#'   argument in either [updateMarkers()] or [deconvolute()]. The elements of
#'   each vector are the values to be tuned for each parameter.
#' @param output Character value, either `"output"` or `"percent"` specifying
#'   which output from the subclass results element resulting from a call to
#'   [deconvolute()]. This deconvolution result is compared against the actual
#'   sample cell numbers in `samples`, using [Rsq_set()].
#' @param force_intercept Logical whether to force intercept through 0.
#' @param verbose Logical whether to show progress.
#' @param ... Optional arguments passed to [deconvolute()] to control fixed
#'   settings.
#' @returns Dataframe with class `'tune_deconv'` whose columns include: the
#'   parameters being tuned via `grid`, cell subclass and R squared.
#' @seealso [plot_tune()] [summary.tune_deconv()]
#' @importFrom stats aggregate
#' @export
tune_deconv <- function(cm, test, samples, grid,
                        output = "output",
                        force_intercept = FALSE,
                        verbose = TRUE, ...) {
  params <- names(grid)
  arg_set1 <- names(formals(updateMarkers))
  arg_set2 <- names(formals(deconvolute))
  if (any(!params %in% c(arg_set1, arg_set2)))
    stop("unknown tuning parameter in `grid`")
  w1 <- which(params %in% arg_set1)
  w2 <- which(params %in% arg_set2)
  grid2 <- if (length(w2) > 0) expand.grid(grid[w2]) else NULL
  if (verbose) {
    message("Tuning parameters: ", paste(params, collapse = ", "))
  }
  
  if (length(w1) > 0) {
    grid1 <- expand.grid(grid[w1])
    if (verbose) pb <- txtProgressBar2()
    res <- lapply(seq_len(nrow(grid1)), function(i) {
      if (verbose) setTxtProgressBar(pb, i / nrow(grid1))
      args <- list(object = cm)
      grid1_row <- grid1[i, , drop = FALSE]
      args <- c(args, grid1_row)
      cm_update <- do.call("updateMarkers", args) |> suppressMessages()
      df2 <- tune_dec(cm_update, test, samples, grid2, output, force_intercept,
                      ...)
      data.frame(grid1_row, df2, row.names = NULL)
    })
    res <- do.call(rbind, res)
    if (verbose) close(pb)
  } else {
    # null grid1
    if (is.null(grid2)) stop("No parameters to tune")
    res <- tune_dec(cm, test, samples, grid2, output, force_intercept, ...)
  }
  
  mres <- aggregate(res$Rsq, by = res[, params, drop = FALSE], FUN = mean,
                    na.rm = TRUE)
  colnames(mres)[which(colnames(mres) == "x")] <- "mean.Rsq"
  w <- which.max(mres$mean.Rsq)
  best_tune <- mres[w, ]
  if (verbose) {
    cat("Best tune:\n")
    print(best_tune, row.names = FALSE, digits = max(3, getOption("digits")-3),
          print.gap = 2L)
  }
  class(res) <- c("tune_deconv", class(res))
  res
}

# tune inner grid of arguments for deconvolute()
tune_dec <- function(cm, test, samples, grid2, output, force_intercept, ...) {
  if (is.null(grid2)) {
    fit <- deconvolute(cm, test, ...) |> suppressMessages()
    fit_output <- fit$subclass[[output]]
    out <- Rsq_set(samples, fit_output, force_intercept)
    df <- data.frame(subclass = names(out), Rsq = out, row.names = NULL)
    return(df)
  }
  # loop grid2
  res <- lapply(seq_len(nrow(grid2)), function(i) {
    dots <- list(...)
    grid2_row <- grid2[i, , drop = FALSE]
    args <- list(mk = cm, test = test)
    args <- c(args, grid2_row)
    if (length(dots)) args[names(dots)] <- dots
    fit <- do.call("deconvolute", args) |> suppressMessages()
    fit_output <- fit$subclass[[output]]
    out <- Rsq_set(samples, fit_output, force_intercept)
    df <- data.frame(grid2_row, subclass = names(out), Rsq = out,
                     row.names = NULL)
  })
  do.call(rbind, res)
}


#' Summarising deconvolution tuning
#' 
#' `summary` method for class `'tune_deconv'`.
#' 
#' @param object dataframe of class `'tune_deconv'`.
#' @param ... further arguments passed to other methods.
#' @returns Prints the row representing the best tuning of parameters (maximum
#'   mean R squared, averaged across subclasses). Invisibly returns a dataframe
#'   of mean R squared values averaged over subclasses.
#' @export
summary.tune_deconv <- function(object, ...) {
  params <- colnames(object)
  params <- params[!params %in% c("subclass", "Rsq")]
  mres <- aggregate(object$Rsq, by = object[, params, drop = FALSE], FUN = mean,
                    na.rm = TRUE)
  colnames(mres)[which(colnames(mres) == "x")] <- "mean.Rsq"
  w <- which.max(mres$mean.Rsq)
  best_tune <- mres[w, ]
  cat("Best tune:\n")
  print(best_tune, row.names = FALSE, digits = max(3, getOption("digits")-3),
        print.gap = 2L)
  invisible(mres)
}


#' Plot tuning curves
#' 
#' Produces a ggplot2 plot of R-squared values generated by [tune_deconv()].
#' 
#' @param result Dataframe of tuning results generated by [tune_deconv()].
#' @param group Character value specifying column in `result` to be grouped by
#'   colour; or `NULL` to average tuning values across the grid and show the
#'   mean effect of varying the parameter specified by `xvar`.
#' @param xvar Character value specifying column in `result` to vary along the x
#'   axis.
#' @param title Character value for the plot title.
#' @returns ggplot2 scatter plot.
#' @details
#' If `group` is `NULL`, the tuning parameter specified by `xvar` is varied on
#' the x axis and R-squared values are averaged over the whole grid to give the
#' mean effect of varying the `xvar` parameter.
#' 
#' If `group` is set to `"subclass"`, then the tuning parameter specified by
#' `xvar` is varied on the x axis. Any other tuning parameters (i.e. if 2 or
#' more have been tuned) are fixed to their best tuned values.
#' 
#' If `group` is set to a different column than `"subclass"`, then the mean
#' R-squared values in `result` are averaged over subclasses. This makes it
#' easier to compare the overall effect (mean R-squared) of 2 tuned parameters
#' which are specified by `xvar` and `group`. Any remaining parameters not shown
#' are fixed to their best tuned values.
#' 
#' @importFrom dplyr near
#' @importFrom ggplot2 geom_line ggtitle mean_se stat_summary theme_bw
#' @export
plot_tune <- function(result, group = "subclass", xvar = colnames(result)[1],
                      title = NULL) {
  params <- colnames(result)
  params <- params[!params %in% c("subclass", "Rsq")]
  if (!xvar %in% params) stop("incorrect `xvar`")
  
  if (is.null(group)) {
    xdiff <- diff(range(result[, xvar], na.rm = TRUE))
    
    p <- ggplot(result, aes(x = .data[[xvar]], y = .data$Rsq)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", col = "black",
                   width = 0.02 * xdiff) +
      stat_summary(fun = mean, geom = "point", col = "black") +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(size = 9),
            axis.text = element_text(colour = "black"))
    return(p)
  }
  if (!group %in% colnames(result)) stop("incorrect `group`")
  by_params <- c(group, xvar)
  fix_params <- params[!params %in% by_params]
  mres <- aggregate(result$Rsq, by = result[, params, drop = FALSE], FUN = mean,
                    na.rm = TRUE)
  colnames(mres)[which(colnames(mres) == "x")] <- "mean.Rsq"
  w <- which.max(mres$mean.Rsq)
  best_tune <- mres[w, ]
  
  if (group == "subclass") {
    # usual plot
    if (length(fix_params)) {
      # 2 or more params tuned, fix using best_tune
      fix <- lapply(fix_params, function(i) {
        near(result[, i], best_tune[, i])
      })
      p <- paste(paste(fix_params, best_tune[, fix_params], sep = " = "),
                 collapse = ", ")
      message("Fix ", p)
      if (is.null(title)) title <- p
      fix <- do.call(cbind, fix)
      if (ncol(fix) > 1) fix <- rowSums(fix) == ncol(fix)
      result <- result[fix, ]
    }
    xdiff <- diff(range(result[, xvar], na.rm = TRUE))
    
    ggplot(result, aes(x = .data[[xvar]], y = .data$Rsq,
                       color = .data[[group]])) +
      geom_line() +
      geom_point() +
      stat_summary(fun.data = mean_se, geom = "errorbar", col = "black",
                   width = 0.02 * xdiff) +
      stat_summary(fun = mean, geom = "point", col = "black") +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(size = 9),
            axis.text = element_text(colour = "black"))
  } else {
    # mean Rsq over subclasses
    if (length(fix_params)) {
      # 3 or more params tuned, fix using best_tune
      fix <- lapply(fix_params, function(i) {
        near(mres[, i], best_tune[, i])
      })
      p <- paste(paste(fix_params, best_tune[, fix_params], sep = " = "),
                 collapse = ", ")
      message("Fix ", p)
      if (is.null(title)) title <- p
      fix <- do.call(cbind, fix)
      if (ncol(fix) > 1) fix <- rowSums(fix) == ncol(fix)
      mres <- mres[fix, ]
    }
    mres[, group] <- factor(mres[, group])
    ggplot(mres, aes(x = .data[[xvar]], y = .data$mean.Rsq,
                       color = .data[[group]])) +
      geom_line() +
      geom_point() +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(size = 9),
            axis.text = element_text(colour = "black"))
  }
}
