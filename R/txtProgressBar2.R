
# Use `utils::setTxtProgressBar()` to set the progress bar and `close()` to
# close it.
# 
#' @importFrom utils flush.console setTxtProgressBar

txtProgressBar2 <- function(min = 0, max = 1, initial = 0, char = "=",
                            width = NA, title = "") {
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (nw == 0) 
    stop("'char' must have a non-zero width")
  if (is.na(width)) {
    width <- getOption("width")
    nt <- nchar(title)
    if (length(nt) == 0) nt <- 0
    width <- width - 22L - nt
    if (nw > 1) 
      width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have 'max' > 'min'")
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste0("\r", title, "  |", strrep(" ", width + 6)))
    cat(paste(c("\r", title, "  |", rep.int(char, nb),
                rep.int(" ", nw * (width - nb)),
                sprintf("| %3d%%", pc)), collapse = ""))
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    end <- Sys.time()
    cat(paste0("  (", format(end - .start, digits = 3), ")\n"))
    flush.console()
    .killed <<- TRUE
  }
  up <- up3
  up(initial)
  .start <- Sys.time()
  structure(list(getVal = getVal, up = up, kill = kill),
            class = "txtProgressBar")
}


mcProgressBar <- function(value, title = "") {
  width <- getOption("width") - 22L - nchar(title)
  nb <- round(width * value)
  pc <- round(100 * value)
  # standard
  p <- paste(c(title, "  |", rep.int("=", nb), rep.int(" ", width - nb),
               sprintf("| %3d%%", pc)), collapse = "")
  if (Sys.getenv("RSTUDIO") == "1") {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::getThemeInfo()$dark) {
      # colour
      p <- paste(c(title, "  |\\x1b[36m", rep.int("=", nb),
                   rep.int(" ", width - nb),
                   sprintf("\\x1b[37m| %3d%%", pc)), collapse = "")
    }
  }
  over_parallel(p)
}

closeProgress <- function(start, title = "") {
  end <- Sys.time()
  mcProgressBar(1, title)
  message_parallel("  (", format(end - start, digits = 3), ")")
}

# prints using shell echo from inside mclapply when run in Rstudio
cat_parallel <- function(...) {
  if (Sys.getenv("RSTUDIO") != "1") return()
  system(sprintf('echo "%s', paste0(..., '\\c"', collapse = "")))
}

message_parallel <- function(...) {
  if (Sys.getenv("RSTUDIO") != "1") return()
  system(sprintf('echo "%s"', paste0(..., collapse = "")))
}

over_parallel <- function(...) {
  if (Sys.getenv("RSTUDIO") != "1") return()
  p <- paste0('\\r', ..., '\\c"', collapse = "")
  system(sprintf('echo "%s', p))
}


