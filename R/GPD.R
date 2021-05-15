#' Create a GPD distribution (Generalized Pareto Distribution)
#'
#'
#' @param loc The location parameter, often noted \eqn{\mu}, which represents
#'   the threshold above which data are taken. Can be any real number. Defaults to `0`.
#' @param scale The scale parameter, often noted \eqn{\sigma}. Can be any
#'   positive number. Defaults to `1`.
#' @param shape The shape parameter, often noted \eqn{\xi}. Can be any
#'   real number. Defaults to `0`, in which case the GPD distribution
#'   correspond to the Exponential distribution. Note that the 'statistical'
#'   convention (also used in package evd, in Wikipedia and in many textbooks)
#'   is used here: \eqn{\xi > 0} correspond to heavy-tailed distributions.
#'   \eqn{\xi < 0} correspond to light-tailed, right-bounded distributions.
#'   Be aware that the opposite convention is used in some disciplines (typically, Hydrology).
#'
#' @return A `GPD` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GPD(loc=0,scale=1,shape=0.2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 0.7)
#' log_pdf(X, 0.7)
#'
#' cdf(X, 0.7)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
GPD <- function(loc = 0, scale = 1, shape = 0) {
  d <- list(loc = loc, scale = scale, shape = shape)
  class(d) <- c("GPD", "distribution")
  d
}

#' @export
print.GPD <- function(x, ...) {
  cat(glue::glue("GPD distribution (loc = {x$loc}, scale = {x$scale}, shape = {x$shape})\n"))
}

#' Draw a random sample from a GPD distribution
#'
#' @inherit GPD examples
#'
#' @param d A `GPD` object created by a call to [GPD()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.GPD <- function(d, n = 1L, ...) {
  if(!feas.GPD(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    out <- evd::rgpd(n,d$loc,d$scale,d$shape)
  }
  out
}

#' Evaluate the probability density function of a GEV distribution
#'
#' @inherit GPD examples
#' @inheritParams random.GPD
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.GPD <- function(d, x, ...) {
  if(!feas.GPD(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::dgpd(x,d$loc,d$scale,d$shape)
  }
  out
}

#' @rdname pdf.GPD
#' @export
#'
log_pdf.GPD <- function(d, x, ...) {
  log(pdf.GPD(d,x))
}

#' Evaluate the cumulative distribution function of a GEV distribution
#'
#' @inherit GPD examples
#' @inheritParams random.GPD
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GPD <- function(d, x, ...) {
  if(!feas.GPD(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::pgpd(x,d$loc,d$scale,d$shape)
  }
  out
}

#' Determine quantiles of a GPD distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GPD examples
#' @inheritParams random.GPD
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.GPD <- function(d, p, ...) {
  if(!feas.GPD(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out=p
    out[p<0] <- NaN
    out[p>1] <- NaN
    out[p==1] <- ifelse(d$shape>=0,Inf,d$loc-d$scale/d$shape)
    out[p==0] <- d$loc
    mask <- (p>0 & p<1)
    out[mask] <- evd::qgpd(p[mask],d$loc,d$scale,d$shape)
  }
  out
}

#' @export
feas.GPD<-function(d){
  d$scale>0
}
