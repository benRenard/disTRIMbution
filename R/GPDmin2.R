#' Create a GPDmin2 distribution (Generalized Pareto Distribution for minima with alternative parameterization)
#'
#'
#' @param loc The location parameter, often noted \eqn{\mu}, which represents
#'   the threshold below which data are taken. Can be any real number.
#'   Defaults to `0`.
#' @param scale The scale parameter, often noted \eqn{\sigma}. Can be any
#'   positive number. Defaults to `1`.
#' @param lbound The lower bound parameter \eqn{\alpha}. Can be any
#'   real number. Defaults to `-1`. It should be strictly smaller than the location
#'   parameter \eqn{\mu}. Relation with the standard GPD parameterization is
#'   \eqn{\alpha=\mu+\sigma/\xi} and \eqn{\xi<0}.
#'
#' @return A `GPDmin2` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GPDmin2(loc=2,scale=1,lbound=0)
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
GPDmin2 <- function(loc = 0, scale = 1, lbound = -1) {
  d <- list(loc = loc, scale = scale, lbound = lbound)
  class(d) <- c("GPDmin2", "distribution")
  d
}

#' @export
print.GPDmin2 <- function(x, ...) {
  cat(glue::glue("GPDmin2 distribution (loc = {x$loc}, scale = {x$scale}, lbound = {x$lbound})\n"))
}

#' Draw a random sample from a GPDmin2 distribution
#'
#' @inherit GPDmin2 examples
#'
#' @param d A `GPDmin2` object created by a call to [GPDmin2()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.GPDmin2 <- function(d, n = 1L, ...) {
  if(!feas.GPDmin2(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    out <- -1*evd::rgpd(n,-1*d$loc,d$scale,getShape(d))
  }
  out
}

#' Evaluate the probability density function of a GPDmin2 distribution
#'
#' @inherit GPDmin2 examples
#' @inheritParams random.GPDmin2
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.GPDmin2 <- function(d, x, ...) {
  if(!feas.GPDmin2(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- evd::dgpd(-1*x,-1*d$loc,d$scale,getShape(d))
  }
  out
}

#' @rdname pdf.GPDmin2
#' @export
#'
log_pdf.GPDmin2 <- function(d, x, ...) {
  log(pdf.GPDmin2(d,x))
}

#' Evaluate the cumulative distribution function of a GPDmin2 distribution
#'
#' @inherit GPDmin2 examples
#' @inheritParams random.GPDmin2
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GPDmin2 <- function(d, x, ...) {
  if(!feas.GPDmin2(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- 1-evd::pgpd(-1*x,-1*d$loc,d$scale,getShape(d))
  }
  out
}

#' Determine quantiles of a GPDmin2 distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GPDmin2 examples
#' @inheritParams random.GPDmin2
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.GPDmin2 <- function(d, p, ...) {
  if(!feas.GPDmin2(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out=p
    out[p<0] <- NaN
    out[p>1] <- NaN
    out[p==1] <- d$loc
    out[p==0] <- d$lbound
    mask <- (p>0 & p<1)
    out[mask] <- -1*evd::qgpd(1-p[mask],-1*d$loc,d$scale,getShape(d))
  }
  out
}

#' @export
feas.GPDmin2<-function(d){
  d$scale>0 & d$lbound<d$loc
}

