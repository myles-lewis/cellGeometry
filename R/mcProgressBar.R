
## Example
# cores <- 4
# my_fun <- function() {
#   start <- Sys.time()
#   out <- mclapply(1:20, function(x) {
#     mcProgressBar(x, 20, cores, "Calculating")
#     Sys.sleep(0.5 + runif(1))
#     rnorm(1e4)
#   }, mc.cores = cores)
#   closeProgress(start, "Calculating")
#   out
# }
# x <- my_fun()

mcProgressBar <- function(val, len = 1L, cores = 1L, title = "",
                          spinner = TRUE) {
  if (val %% cores != 0) return(if (spinner) mcSpinner(val, title))
  width <- getOption("width") - 22L - nchar(title)
  nb <- round(width * val / len)
  pc <- round(100 * val / len)
  i <- val %% 4 +1
  sp <- if (pc == 100 | !spinner) "  " else c("/ ", "- ", "\\\ ", "| ")[i]
  if (title != "") title <- paste0(title, " ")
  # standard
  p <- paste(c(title, sp, "|", rep.int("=", nb), rep.int(" ", width - nb),
               sprintf("| %3d%%", pc)), collapse = "")
  if (Sys.getenv("RSTUDIO") == "1") {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::getThemeInfo()$dark) {
      # colour
      p <- paste(c("\\x1b[37m", title, sp, "|\\x1b[36m", rep.int("=", nb),
                   rep.int(" ", width - nb),
                   sprintf("\\x1b[37m| %3d%%", pc)), collapse = "")
    }
  }
  over_parallel(p)
}

mcSpinner <- function(val, title) {
  i <- val %% 4 +1
  sp <- c("/", "-", "\\\ ", "|")[i]
  if (title != "") title <- paste0(title, " ")
  over_parallel(title, sp)
}

closeProgress <- function(start, title = "") {
  end <- Sys.time()
  mcProgressBar(1, title = title)
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
