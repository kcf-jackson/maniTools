#' A function to run the shiny app
#' @export
run_example <- function() {
  appDir <- system.file("shiny-examples", "maniTools", package = "maniTools")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `maniTools`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
