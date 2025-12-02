#' SHINY MLE
#'
#' @returns shiny app
#' @export
#' @import shiny
#'
#' @examples
#' \dontrun{shinymle()}
shinymle = function(){
  shiny::runApp(system.file("SHINY", package = "MATH4753Schmidt2025"), launch.browser = TRUE)
}
