#' Create a GPDmin distribution (Generalized Pareto Distribution for minima)
#'
#'
#' @param loc The location parameter, often noted \eqn{\mu}, which represents
#'   the threshold below which data are taken. Can be any real number. Defaults to `0`.
#' @param scale The scale parameter, often noted \eqn{\sigma}. Can be any
#'   positive number. Defaults to `1`.
#' @param shape The shape parameter, often noted \eqn{\xi}. Can be any
#'   real number. Defaults to `0`. Note that the convention used here derives from
#'   the one used for the GPD distribution:
#'   \eqn{\xi > 0} leads to a right-bounded distribution.
#'   \eqn{\xi < 0} leads to a left-bounded distribution.
#'
#' @return A `GPDmin` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GPDmin(loc=1,scale=1,shape=-0.2)
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
GPDmin <- function(loc = 0, scale = 1, shape = 0) {
  d <- list(loc = loc, scale = scale, shape = shape)
  class(d) <- c("GPDmin", "distribution")
  d
}

#' @export
print.GPDmin <- function(x, ...) {
  cat(glue::glue("GPDmin distribution (loc = {x$loc}, scale = {x$scale}, shape = {x$shape})\n"))
}

#' Draw a random sample from a GPDmin distribution
#'
#' @inherit GPDmin examples
#'
#' @param d A `GPDmin` object created by a call to [GPDmin()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.GPDmin <- function(d, n = 1L, ...) {
  if(!feas.GPDmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    out <- -1*evd::rgpd(n,-1*d$loc,d$scale,d$shape)
  }
  out
}

#' Evaluate the probability density function of a GPDmin distribution
#'
#' @inherit GPDmin examples
#' @inheritParams random.GPDmin
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.GPDmin <- function(d, x, ...) {
  if(!feas.GPDmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::dgpd(-1*x,-1*d$loc,d$scale,d$shape)
  }
  out
}

#' @rdname pdf.GPDmin
#' @export
#'
log_pdf.GPDmin <- function(d, x, ...) {
  log(pdf.GPDmin(d,x))
}

#' Evaluate the cumulative distribution function of a GPDmin distribution
#'
#' @inherit GPDmin examples
#' @inheritParams random.GPDmin
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GPDmin <- function(d, x, ...) {
  if(!feas.GPDmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- 1-evd::pgpd(-1*x,-1*d$loc,d$scale,d$shape)
  }
  out
}

#' Determine quantiles of a GPDmin distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GPDmin examples
#' @inheritParams random.GPDmin
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.GPDmin <- function(d, p, ...) {
  if(!feas.GPDmin(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out=p
    out[p<0] <- NaN
    out[p>1] <- NaN
    out[p==1] <- d$loc
    out[p==0] <- ifelse(d$shape>=0,-Inf,d$loc+d$scale/d$shape)
    mask <- (p>0 & p<1)
    out[mask] <- -1*evd::qgpd(1-p[mask],-1*d$loc,d$scale,d$shape)
  }
  out
}

#' @export
feas.GPDmin<-function(d){
  d$scale>0
}
