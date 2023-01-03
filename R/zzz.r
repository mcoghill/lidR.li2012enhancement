#' @useDynLib lidR_li2012_enhancement
#' @importFrom Rcpp sourceCpp
#' @import data.table
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("lidR_li2012_enhancement", libpath)
}
