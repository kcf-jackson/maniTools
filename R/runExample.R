#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "maniTools", package = "maniTools")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `maniTools`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
