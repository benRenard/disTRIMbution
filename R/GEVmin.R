#' Create a GEVmin distribution (Generalized Extreme Value for minima)
#'
#'
#' @param loc The location parameter, often noted \eqn{\mu}. Can be any real number.
#'   Defaults to `0`.
#' @param scale The scale parameter, often noted \eqn{\sigma}. Can be any
#'   positive number. Defaults to `1`.
#' @param shape The shape parameter, often noted \eqn{\xi}. Can be any
#'   real number. Defaults to `0`. Note that the convention used here derives from
#'   the one used for the GEV distribution:
#'   \eqn{\xi > 0} leads to a right-bounded distribution.
#'   \eqn{\xi < 0} leads to a left-bounded distribution.
#'
#' @return A `GEVmin` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GEVmin(loc=0,scale=1,shape=-0.2)
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
GEVmin <- function(loc = 0, scale = 1, shape = 0) {
  d <- list(loc = loc, scale = scale, shape = shape)
  class(d) <- c("GEVmin", "distribution")
  d
}

#' @export
print.GEVmin <- function(x, ...) {
  cat(glue::glue("GEVmin distribution (loc = {x$loc}, scale = {x$scale}, shape = {x$shape})\n"))
}

#' Draw a random sample from a GEVmin distribution
#'
#' @inherit GEVmin examples
#'
#' @param d A `GEVmin` object created by a call to [GEVmin()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.GEVmin <- function(d, n = 1L, ...) {
  if(!feas.GEVmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    out <- -1*evd::rgev(n,-1*d$loc,d$scale,d$shape)
  }
  out
}

#' Evaluate the probability density function of a GEVmin distribution
#'
#' @inherit GEVmin examples
#' @inheritParams random.GEVmin
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.GEVmin <- function(d, x, ...) {
  if(!feas.GEVmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::dgev(-1*x,-1*d$loc,d$scale,d$shape)
  }
  out
}

#' @rdname pdf.GEVmin
#' @export
#'
log_pdf.GEVmin <- function(d, x, ...) {
  log(pdf.GEVmin(d,x))
}

#' Evaluate the cumulative distribution function of a GEVmin distribution
#'
#' @inherit GEVmin examples
#' @inheritParams random.GEVmin
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GEVmin <- function(d, x, ...) {
  if(!feas.GEVmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- 1-evd::pgev(-1*x,-1*d$loc,d$scale,d$shape)
  }
  out
}

#' Determine quantiles of a GEVmin distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GEVmin examples
#' @inheritParams random.GEVmin
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.GEVmin <- function(d, p, ...) {
  if(!feas.GEVmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out=p
    out[p<0] <- NaN
    out[p>1] <- NaN
    out[p==1] <- ifelse(d$shape<=0,Inf,d$loc+d$scale/d$shape)
    out[p==0] <- ifelse(d$shape>=0,-Inf,d$loc+d$scale/d$shape)
    mask <- (p>0 & p<1)
    out[mask] <- -1*evd::qgev(1-p[mask],-1*d$loc,d$scale,d$shape)
  }
  out
}

#' @export
feas.GEVmin<-function(d){
  d$scale>0
}
