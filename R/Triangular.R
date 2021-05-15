#' Create a Triangular distribution
#'
#'
#' @param a The a parameter (lower bound). `a` can be any value in the set of real
#'   numbers. Defaults to `0`.
#' @param b The b parameter (upper bound). `b` can be any value in the set of real
#'   numbers. It should be strictly bigger than `a`. Defaults to `1`.
#' @param c The c parameter (peak). `c` can be any value in the set of real
#'   numbers. It should be strictly between `a` and `b`. Defaults to `0.5`.
#'
#' @return A `Triangular` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Triangular(a=0,b=4,c=1)
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
Triangular <- function(a = 0, b = 1, c = 0.5) {
  d <- list(a = a, b = b, c = c)
  class(d) <- c("Triangular", "distribution")
  d
}

#' @export
print.Triangular <- function(x, ...) {
  cat(glue::glue("Triangular distribution (a = {x$a}, b = {x$b}, c = {x$c})\n"))
}

#' Draw a random sample from a Triangular distribution
#'
#' @inherit Triangular examples
#'
#' @param d A `Triangular` object created by a call to [Triangular()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Triangular <- function(d, n = 1L, ...) {
  if(!feas.Triangular(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    u <- stats::runif(n)
    out <- quantile.Triangular(d,u)
  }
  out
}

#' Evaluate the probability density function of a Triangular distribution
#'
#' @inherit Triangular examples
#' @inheritParams random.Triangular
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.Triangular <- function(d, x, ...) {
  if(!feas.Triangular(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- x
    mask <- x<d$a; out[mask] <- 0
    mask <- (x>=d$a & x<d$c); out[mask] <- (2*(x[mask]-d$a))/((d$b-d$a)*(d$c-d$a))
    mask <- (x>=d$c & x<d$b); out[mask] <- (2*(d$b-x[mask]))/((d$b-d$a)*(d$b-d$c))
    mask <- x>=d$b; out[mask] <- 0
  }
  out
}

#' @rdname pdf.Triangular
#' @export
#'
log_pdf.Triangular <- function(d, x, ...) {
  log(pdf.Triangular(d,x))
}

#' Evaluate the cumulative distribution function of a Triangular distribution
#'
#' @inherit Triangular examples
#' @inheritParams random.Triangular
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Triangular <- function(d, x, ...) {
  if(!feas.Triangular(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- x
    mask <- x<d$a; out[mask] <- 0
    mask <- (x>=d$a & x<d$c); out[mask] <- ((x[mask]-d$a)^2)/((d$b-d$a)*(d$c-d$a))
    mask <- (x>=d$c & x<d$b); out[mask] <- 1-((d$b-x[mask])^2)/((d$b-d$a)*(d$b-d$c))
    mask <- x>=d$b; out[mask] <- 1
  }
  out
}

#' Determine quantiles of a Triangular distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Triangular examples
#' @inheritParams random.Triangular
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Triangular <- function(d, p, ...) {
  if(!feas.Triangular(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out <- p
    out[p<0] <- NaN
    out[p>1] <- NaN
    plim <- (d$c-d$a)/(d$b-d$a)
    mask <- (p>=0 & p<plim);out[mask] <- d$a+sqrt(p[mask]*(d$b-d$a)*(d$c-d$a))
    mask <- (p>=plim & p<=1);out[mask] <- d$b-sqrt((1-p[mask])*(d$b-d$a)*(d$b-d$c))
  }
  out
}

#' @export
feas.Triangular<-function(d){
  if( d$b<=d$a | d$c<d$a | d$c>d$b ){out=FALSE} else {out=TRUE}
  return(out)
}
