#' Calculate geometric mean
#'
#' Credit to:
#' http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
#'
#' @param x Numeric vector
#' @param na.rm Ignore NAs or not
#'
#' @return Single value
#' @export
#'
#' @examples
#' gmean(c(1, 2, 3, 4, 5, NA))
gmean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
