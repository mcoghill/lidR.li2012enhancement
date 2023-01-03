#' @useDynLib lidR.li2012enhancement
#' @importFrom Rcpp sourceCpp
#' @import data.table
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("lidR.li2012enhancement", libpath)
}
