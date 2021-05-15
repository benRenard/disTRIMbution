#***************************************************************************----
# Generics ----

# Import generics used in disTRIMbution - avoids R CMD check warnings / notes
#' @importFrom distributions3 pdf log_pdf cdf random quantile
foo <- NULL # useless - just a hack to avoid roxygen misinterpretation

#' Check feasability of a distribution (e.g. positive scale,
#' lower / upper bounds in correct order, etc.)
#'
#' @param d distributions3 object
#' @return a logical, is the distribution feasible?
#' @examples
#' feas(Triangular())
#' feas(Triangular(a=2,b=1))
#' @export
feas <- function(d) {
  UseMethod("feas")
}


#***************************************************************************----
# Misc. ----

#' Get left / middle / right probabilities
#'
#' @param d distributions3 object
#' @param trim Numeric vector of size 2, trimming bounds.
#' @return a list with values left,right,mid.
#' @keywords internal

getPs <- function(d,trim=c(-Inf,Inf)){
  pleft=cdf_trim(d,trim[1])
  pright=1-cdf_trim(d,trim[2])
  pmid=1-pleft-pright
  if(!is.na(pmid)){if(pmid<0) pmid=NA}
  return(list(left=pleft,right=pright,mid=pmid))
}
