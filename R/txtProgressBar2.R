
# Use `utils::setTxtProgressBar()` to set the progress bar and `close()` to
# close it.
# 
#' @importFrom utils flush.console setTxtProgressBar

txtProgressBar2 <- function(min = 0, max = 1, initial = 0, char = "=",
                            width = NA, title = "", eta = TRUE) {
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
    if (eta) {
      curr <- Sys.time()
      dur <- as.numeric(difftime(curr, .start, units = "secs"))
      rem <- format_dur((100 - pc) / pc * dur)
      tim <- if (pc != 100) paste("  eta", rem) else ""
      tim <- str_pad(tim, 15)
    }
    cat(paste0("\r", title, "  |", strrep(" ", width + 6)))
    cat(paste(c("\r", title, "  |", strrep(char, nb),
                strrep(" ", nw * (width - nb)),
                sprintf("| %3d%%", pc), tim), collapse = ""))
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


format_dur <- function(x) {
  x <- as.numeric(x)
  s <- x %% 60
  m <- floor(x / 60)
  h <- floor(m / 60)
  if (h > 0) {
    m <- m %% 60
    return(paste0(h, "h ", m, "m"))
  }
  if (m > 0) {
    return(paste0(m, "m ", floor(s), "s"))
  }
  s <- format(round(s, 1), digits = 2)
  paste0(s, " secs")
}


str_pad <- function(x, n) {
  nc <- nchar(x)
  if (nc < n) x <- paste0(x, strrep(" ", n - nc))
  x
}
