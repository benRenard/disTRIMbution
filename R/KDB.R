#' Create a KDB distribution (Kuramaswamy double-bounded)
#'
#'
#' @param alpha The alpha parameter (first shape parameter).
#'   `alpha` can be any value strictly greater than zero. Defaults to `1`.
#' @param beta The beta parameter (second shape parameter).
#'   `beta` can be any value strictly greater than zero. Defaults to `1`.
#'
#' @return A `KDB` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- KDB(alpha=2,beta=5)
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
KDB <- function(alpha = 1, beta = 1) {
  d <- list(alpha = alpha, beta = beta)
  class(d) <- c("KDB", "distribution")
  d
}

#' @export
print.KDB <- function(x, ...) {
  cat(glue::glue("KDB distribution (alpha = {x$alpha}, beta = {x$beta})\n"))
}

#' Draw a random sample from a KDB distribution
#'
#' @inherit KDB examples
#'
#' @param d A `KDB` object created by a call to [KDB()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.KDB <- function(d, n = 1L, ...) {
  if(!feas.KDB(d)){ # unfeasible parameters
    out <- return(rep(NaN,n))
  } else {
    u <- stats::runif(n)
    out <- quantile.KDB(d,u)
  }
  out
}

#' Evaluate the probability density function of a KDB distribution
#'
#' @inherit KDB examples
#' @inheritParams random.KDB
#'
#' @param x A vector of elements whose pdf you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of pdf values, one for each element of `x`.
#' @export
#'
pdf.KDB <- function(d, x, ...) {
  if(!feas.KDB(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- x
    mask <- (x>=0 & x<=1)
    out[!mask] <- 0
    out[mask] <- d$alpha*d$beta * x[mask]^(d$alpha-1) * (1-x[mask]^d$alpha)^(d$beta-1)
  }
  out
}

#' @rdname pdf.KDB
#' @export
#'
log_pdf.KDB <- function(d, x, ...) {
  log(pdf.KDB(d,x))
}

#' Evaluate the cumulative distribution function of a KDB distribution
#'
#' @inherit KDB examples
#' @inheritParams random.KDB
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.KDB <- function(d, x, ...) {
  if(!feas.KDB(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(x)))
  } else {
    out <- x
    out[x<0] <- 0
    out[x>1] <- 1
    mask <- (x>=0 & x<=1)
    out[mask] <- 1-(1-x[mask]^d$alpha)^d$beta
  }
  out
}

#' Determine quantiles of a KDB distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit KDB examples
#' @inheritParams random.KDB
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.KDB <- function(d, p, ...) {
  if(!feas.KDB(d)){ # unfeasible parameters
    out <- return(rep(NaN,length(p)))
  } else {
    out <- p
    out[p<0] <- NaN
    out[p>1] <- NaN
    mask <- (p>=0 & p<=1)
    out[mask] <- (1-(1-p[mask])^(1/d$beta))^(1/d$alpha)
  }
  out
}

#' @export
feas.KDB<-function(d){
  (d$alpha>0) & (d$beta>0)
}
