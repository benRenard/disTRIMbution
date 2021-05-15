#***************************************************************************----
# get/set information from/to distributions3 object ----

#' Get parameter vector
#'
#' Get parameter vector from a distributions3 object.
#' Note that this function is highly dependent to the
#' structure of a distributions3 object. If one day the developpers decide
#' to add more fields to this structure, this function will have to be changed!
#'
#' @param d distributions3 object
#' @return a numeric vector containing the parameters of d.
#' @examples
#' \dontrun{
#' d <- Normal()
#' getParVector(d)
#' }
#' @keywords internal

getParVector <- function(d){
  as.numeric(unclass(d))
}

#' Get parameter names
#'
#' Get parameter names from a distributions3 object.
#'
#' @param d distributions3 object
#' @return a character vector containing the parameter names of d.
#' @examples
#' \dontrun{
#' d <- Normal()
#' getParName(d)
#' }
#' @keywords internal

getParName <- function(d){
  names(unclass(d))
}

#' copyDist
#'
#' Make a copy of a distributions3 object and reset its parameter vector.
#'
#' @param d distributions3 object to be copied (provides the distribution family)
#' @param param numeric vector of size npar(d), new parameter vector
#' @return a distributions3 object.
#' @examples
#' \dontrun{
#' d1 <- Normal()
#' d2 <- copyDist(d,param=c(10,5))
#' }
#' @keywords internal

copyDist <- function(d,param=getParVector(d)){
  out <- d
  out[] <- param
  return(out)
}

#' Get number of parameters
#'
#' @param d distributions3 object
#' @return the number of parameters of d.
#' @examples
#' \dontrun{
#' d <- Normal()
#' nPar(d)
#' }
#' @keywords internal

nPar <- function(d){
  length(getParVector(d))
}
