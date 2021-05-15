#' Create a GEV distribution (Generalized Extreme Value)
#'
#'
#' @param loc The location parameter, often noted \eqn{\mu}. Can be any real number.
#'   Defaults to `0`.
#' @param scale The scale parameter, often noted \eqn{\sigma}. Can be any
#'   positive number. Defaults to `1`.
#' @param shape The shape parameter, often noted \eqn{\xi}. Can be any
#'   real number. Defaults to `0`, in which case the GEV distribution
#'   is also known as the Gumbel distribution. Note that the 'statistical'
#'   convention (also used in package evd, in Wikipedia and in many textbooks)
#'   is used here: \eqn{\xi > 0} correspond to heavy-tailed, left-bounded
#'   distributions (aka Frechet Family or Type-II distributions).
#'   \eqn{\xi < 0} correspond to light-tailed, right-bounded
#'   distributions (aka Weibull Family or Type-III distributions). Be aware that
#'   the opposite convention is used in some disciplines (typically, Hydrology).
#'
#' @return A `GEV` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GEV(loc=0,scale=1,shape=0.2)
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
GEV <- function(loc = 0, scale = 1, shape = 0) {
  d <- list(loc = loc, scale = scale, shape = shape)
  class(d) <- c("GEV", "distribution")
  d
}

#' @export
print.GEV <- function(x, ...) {
  cat(glue::glue("GEV distribution (loc = {x$loc}, scale = {x$scale}, shape = {x$shape})\n"))
}

#' Draw a random sample from a GEV distribution
#'
#' @inherit GEV examples
#'
#' @param d A `GEV` object created by a call to [GEV()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.GEV <- function(d, n = 1L, ...) {
  if(!feas.GEV(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    out <- evd::rgev(n,d$loc,d$scale,d$shape)
  }
  out
}

#' Evaluate the probability density function of a GEV distribution
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.GEV <- function(d, x, ...) {
  if(!feas.GEV(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::dgev(x,d$loc,d$scale,d$shape)
  }
  out
}

#' @rdname pdf.GEV
#' @export
#'
log_pdf.GEV <- function(d, x, ...) {
  log(pdf.GEV(d,x))
}

#' Evaluate the cumulative distribution function of a GEV distribution
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GEV <- function(d, x, ...) {
  if(!feas.GEV(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::pgev(x,d$loc,d$scale,d$shape)
  }
  out
}

#' Determine quantiles of a GEV distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.GEV <- function(d, p, ...) {
  if(!feas.GEV(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out=p
    out[p<0] <- NaN
    out[p>1] <- NaN
    out[p==1] <- ifelse(d$shape>=0,Inf,d$loc-d$scale/d$shape)
    out[p==0] <- ifelse(d$shape<=0,-Inf,d$loc-d$scale/d$shape)
    mask <- (p>0 & p<1)
    out[mask] <- evd::qgev(p[mask],d$loc,d$scale,d$shape)
  }
  out
}

#' @export
feas.GEV<-function(d){
  d$scale>0
}
