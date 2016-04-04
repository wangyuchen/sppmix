#' @useDynLib sppmix
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function(libpath) {
  library.dynam.unload("sppmix", libpath)
}


#' @importFrom spatstat owin
#' @export
spatstat::owin


#' @importFrom spatstat square
#' @export
spatstat::square
