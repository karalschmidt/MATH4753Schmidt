#' Birthday Function
#'
#' @param x Random sample of size x of people
#'
#' @returns The probability of two of more people in x having a common birthday
#' @export
#'
#' @examples
#' \dontrun{birthday(23,365)}
birthday = function(x){
  1 - exp(lchoose(365,x) + lfactorial(x) - x*log(365))
}
