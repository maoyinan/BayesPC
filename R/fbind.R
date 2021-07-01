#' fuction to...
#'
#' @param a factor
#' @param b factor
#'
#' @return factor
#' @export
#'
#' @examples
fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}
