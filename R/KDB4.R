#' Create a KDB4 distribution (4-parameter Kuramaswamy double-bounded)
#'
#'
#' @param alpha The alpha parameter (first shape parameter).
#'   `alpha` can be any value strictly greater than zero. Defaults to `1`.
#' @param beta The beta parameter (second shape parameter).
#'   `beta` can be any value strictly greater than zero. Defaults to `1`.
#' @param a The a parameter (lower bound). `a` can be any value in the set of real
#'   numbers. Defaults to `0`.
#' @param b The b parameter (upper bound). `b` can be any value in the set of real
#'   numbers. It should be strictly bigger than `a`. Defaults to `1`.
#'
#' @return A `KDB4` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- KDB4(alpha=2,beta=5,a=0.5,b=3)
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
KDB4 <- function(alpha = 1, beta = 1, a = 0, b = 1) {
  d <- list(alpha = alpha, beta = beta, a = a, b = b)
  class(d) <- c("KDB4", "distribution")
  d
}

#' @export
print.KDB4 <- function(x, ...) {
  cat(glue::glue("KDB4 distribution (alpha = {x$alpha}, beta = {x$beta}, a = {x$a}, b = {x$b})\n"))
}

#' Draw a random sample from a KDB4 distribution
#'
#' @inherit KDB4 examples
#'
#' @param d A `KDB4` object created by a call to [KDB4()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.KDB4 <- function(d, n = 1L, ...) {
  if(!feas.KDB4(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    out <- d$a+(d$b-d$a)*random.KDB(d,n)
  }
  out
}

#' Evaluate the probability density function of a KDB4 distribution
#'
#' @inherit KDB4 examples
#' @inheritParams random.KDB4
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.KDB4 <- function(d, x, ...) {
  if(!feas.KDB4(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- pdf.KDB(d,(x-d$a)/(d$b-d$a)) * (1/(d$b-d$a))
  }
  out
}

#' @rdname pdf.KDB4
#' @export
#'
log_pdf.KDB4 <- function(d, x, ...) {
  log(pdf.KDB4(d,x))
}

#' Evaluate the cumulative distribution function of a KDB4 distribution
#'
#' @inherit KDB4 examples
#' @inheritParams random.KDB4
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.KDB4 <- function(d, x, ...) {
  if(!feas.KDB4(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- cdf.KDB(d,(x-d$a)/(d$b-d$a))
  }
  out
}

#' Determine quantiles of a KDB4 distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit KDB4 examples
#' @inheritParams random.KDB4
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.KDB4 <- function(d, p, ...) {
  if(!feas.KDB4(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out <- d$a+(d$b-d$a)*quantile.KDB(d,p)
  }
  out
}

#' @export
feas.KDB4<-function(d){
  (d$alpha>0) & (d$beta>0) & (d$b>d$a)
}
