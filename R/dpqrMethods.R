#***************************************************************************----
# Main exported functions ----

#' Pdf of a trimmed distribution
#'
#' Evaluate the pdf of a trimmed distribution.
#' Note that in the case of a rectified distribution, the
#' pdf at a bound is computed from the parent distribution d,
#' and is hence not equal to the probability of being equal to the bound.
#'
#' @param d distributions3 object
#' @param x Numeric vector, values at which the pdf is evaluated
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @return a numeric vector containing the pdf values.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' # Define a grid of values x
#' x=seq(-1,2,,100)
#' # Compare pdf of basic / rectified / truncated Gaussian distributions
#' y0=pdf_trim(norm,x)
#' yr=pdf_trim(norm,x,trim=c(0,1))
#' yt=pdf_trim(norm,x,trim=c(0,1),doTrunc=TRUE)
#' plot(x,y0,type='l',ylim=c(0,max(yt)))
#' points(x,yr,col='red')
#' points(x,yt,col='blue',pch=3)
#' @export

pdf_trim <- function(d,x,trim=c(-Inf,Inf),doTrunc=FALSE){
  if(trim[2]<=trim[1]){
    message('Fatal: trimming bounds are in decreasing order')
    return(NA)
  }
  out <- distributions3::pdf(d,x)
  if(doTrunc){
    p=getPs(d,trim)
    if(p$mid<0) {out <- NA*out} else {out <- out/p$mid}
  }
  out[x<trim[1] | x>trim[2]] <- 0
  return(out)
}

#' Cdf of a trimmed distribution
#'
#' Evaluate the cdf of a trimmed distribution.
#'
#' @param d distributions3 object
#' @param x Numeric vector, values at which the cdf is evaluated
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @return a numeric vector containing the cdf values.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' # Define a grid of values x
#' x=seq(-1,2,,100)
#' # Compare cdf of basic / rectified / truncated Gaussian distributions
#' y0=cdf_trim(norm,x)
#' yr=cdf_trim(norm,x,trim=c(0,1))
#' yt=cdf_trim(norm,x,trim=c(0,1),doTrunc=TRUE)
#' plot(x,y0,type='l',ylim=c(0,1))
#' points(x,yr,col='red')
#' points(x,yt,col='blue',pch=3)
#' @export

cdf_trim <- function(d,x,trim=c(-Inf,Inf),doTrunc=FALSE){
  if(trim[2]<=trim[1]){
    message('Fatal: trimming bounds are in decreasing order')
    return(NA)
  }
  out <- distributions3::cdf(d,x)
  if(doTrunc){
    p=getPs(d,trim)
    if(p$mid<0) {out <- NA*out} else {out <- (out-p$left)/p$mid}
  }
  out[x<trim[1]] <- 0
  out[x>trim[2]] <- 1
  return(out)
}

#' Quantile of a trimmed distribution
#'
#' Evaluate the quantile of a trimmed distribution.
#'
#' @param d distributions3 object
#' @param p Numeric vector in [0,1], probabilities at which quantiles are evaluated
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @return a numeric vector containing the quantiles.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' # Define a grid of values p
#' p=seq(0,1,,100)
#' # Compare quantiles of basic / rectified / truncated Gaussian distributions
#' y0=quantile_trim(norm,p)
#' yr=quantile_trim(norm,p,trim=c(0,1))
#' yt=quantile_trim(norm,p,trim=c(0,1),doTrunc=TRUE)
#' plot(p,y0,type='l')
#' points(p,yr,col='red')
#' points(p,yt,col='blue',pch=3)
#' @export

quantile_trim <- function(d,p,trim=c(-Inf,Inf),doTrunc=FALSE){
  if(trim[2]<=trim[1]){
    message('Fatal: trimming bounds are in decreasing order')
    return(NA)
  }
  probs=getPs(d,trim)
  if(doTrunc){
    ptrans <- probs$left+p*probs$mid
    out <- distributions3::quantile(d,ptrans)
  } else {
    out <- distributions3::quantile(d,p)
    out[p<probs$left] <- trim[1]
    out[p>1-probs$right] <- trim[2]
  }
  return(out)
}

#' Random numbers from a trimmed distribution
#'
#' Generate random numbers from a trimmed distribution.
#'
#' @param d distributions3 object
#' @param n Integer, number of samples to draw
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @return a numeric vector containing the generated values
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' # Generate values from basic / rectified / truncated Gaussian distributions
#' n=1000
#' y0=random_trim(norm,n)
#' yr=random_trim(norm,n,trim=c(0,1))
#' yt=random_trim(norm,n,trim=c(0,1),doTrunc=TRUE)
#' # plot results
#' par(mfrow=c(2,3))
#' lim=c(-0.75,2.25)
#' plot(y0,ylim=lim);plot(yr,col='red',ylim=lim);plot(yt,col='blue',ylim=lim)
#' hist(y0,xlim=lim);hist(yr,col='red',xlim=lim);hist(yt,col='blue',xlim=lim)
#' @export

random_trim <- function(d,n,trim=c(-Inf,Inf),doTrunc=FALSE){
  if(trim[2]<=trim[1]){
    message('Fatal: trimming bounds are in decreasing order')
    return(NA)
  }
  if(doTrunc){
    p <- stats::runif(n)
    out <- quantile_trim(d,p,trim,TRUE)
  } else {
    out <- distributions3::random(d,n)
    out[out<trim[1]] <- trim[1]
    out[out>trim[2]] <- trim[2]
  }
  return(out)
}
