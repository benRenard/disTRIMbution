#***************************************************************************----
# Main exported functions ----

#' Estimate a trimmed distribution
#'
#' Estimate the parameters of a trimmed distribution using
#' maximum likelihood.
#'
#' @param d distributions3 object, assumed distribution.
#'          The parameters specified in d are interpreted as initial guesses
#'          and should be both feasible (e.g. no negative scale parameter)
#'          and data-compatible (i.e. not lead to a data point having a null pdf)
#' @param y Numeric vector, data
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param par.low Numeric vector of size nPar(d), lower bounds for parameters.
#' @param par.high Numeric vector of size nPar(d), higher bounds for parameters.
#' @param method Character string, optimization method, see ?optim.
#' @param control List, controls for optimization algorithm, see ?optim.
#' @return the estimated distribution as a distributions3 object.
#'     Note that parameters may be NA's if the estimation failed.
#' @examples
#' # EXAMPLE 1
#' # Define Normal distribution
#' norm <- Normal(mu=1,sigma=2)
#' # generate data and reset all negative values to zero
#' y <- random(norm,200); y[y<0] <- 0
#' # Estimate basic / rectified / truncated Gaussian distributions.
#' # The latter returns NA's because y=0 values occur,
#' # which is not compatible with a zero-truncated distribution.
#' estim_ML(d=norm,y=y)
#' estim_ML(d=norm,y=y,trim=c(0,Inf))
#' estim_ML(d=norm,y=y,trim=c(0,Inf),doTrunc=TRUE)
#'
#' # EXAMPLE 2
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.25)
#' # generate data and remove all values outside [0;1]
#' z <- random(norm,200); y <- z[z>0 & z<1]
#' # Estimate basic / rectified / truncated Gaussian distributions.
#' estim_ML(d=norm,y=y)
#' estim_ML(d=norm,y=y,trim=c(0,1))
#' estim_ML(d=norm,y=y,trim=c(0,1),doTrunc=TRUE)
#' # Good practice: explicitly set parameter bounds, here sigma>0
#' estim_ML(d=norm,y=y,par.low=c(-Inf,0))
#' @export

estim_ML <- function(d,y,trim=c(-Inf,Inf),doTrunc=FALSE,
                     par.low=rep(-Inf,nPar(d)),
                     par.high=rep(Inf,nPar(d)),
                     method='Nelder-Mead',
                     control=list(fnscale=-1,maxit=10000)){
  # Check trimming makes sense
  if(trim[2]<=trim[1]){
    message('Fatal: trimming bounds are in decreasing order')
    return(NA)
  }
  # Initialise returned value in case of failure
  fail=copyDist(d,param=rep(NaN,nPar(d)))
  # Check there are no data outside of the trimming bounds
  if(any( y<trim[1] | y>trim[2] )){
    mess=paste('Fatal: all values should be between',trim[1],'and',trim[2],'(included)')
    warning(mess)
    return(fail)
  }
  # Get initial parameter values
  initpar=getParVector(d)
  # optimize
  w=tryCatch(stats::optim(par=initpar,fn=trimLL,gr=NULL,d=d,y=y,
                   trim=trim,doTrunc=doTrunc,
                   par.low=par.low,par.high=par.high,
                   method=method,control=control),
             error=function(e) list(par=fail,convergence=666,value=NA,counts=0)
  )
  if(w$convergence==0){ # worked fine, return parameters
    return(copyDist(d,param=w$par))
  } else { # problem, return NaN's
    return(fail)
  }
}

#' Parametric bootstrap
#'
#' Applies parametric bootstrap to an estimated distribution
#'
#' @param d distributions3 object, estimated distribution.
#'          Typically, the output of estim_ML.
#' @param n Integer, size of each resampled dataset.
#'     Typically length(y), where y are calibration data.
#' @param nsim Integer, number of resampled datasets.
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param par.low Numeric vector of size nPar(d), lower bounds for parameters.
#' @param par.high Numeric vector of size nPar(d), higher bounds for parameters.
#' @param method Character string, optimization method, see ?optim.
#' @param control List, controls for optimization algorithm, see ?optim.
#' @return a matrix of size nsim * nPar(d) containing the estimated parameters
#' @examples
#' # EXAMPLE 1
#' # Define Normal distribution
#' norm <- Normal(mu=1,sigma=2)
#' # generate data and reset all negative values to zero
#' y <- random(norm,200); y[y<0] <- 0
#' # Estimate rectified Gaussian distribution.
#' fit <- estim_ML(d=norm,y=y,trim=c(0,Inf),par.low=c(-Inf,0))
#' # Do parametric bootstrap.
#' boot <- PBoot(fit,length(y),trim=c(0,Inf),par.low=c(-Inf,0),nsim=100)
#' plot(boot)
#'
#' # EXAMPLE 2
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.25)
#' # generate data and remove all values outside [0;1]
#' z <- random(norm,200); y <- z[z>0 & z<1]
#' # Estimate truncated Gaussian distribution.
#' fit <- estim_ML(d=norm,y=y,trim=c(0,1),doTrunc=TRUE,par.low=c(-Inf,0))
#' # Do parametric bootstrap.
#' boot <- PBoot(fit,length(y),trim=c(0,1),doTrunc=TRUE,par.low=c(-Inf,0),nsim=100)
#' plot(boot)
#' @export

PBoot<-function(d,n,nsim=1000,
                trim=c(-Inf,Inf),doTrunc=FALSE,
                par.low=rep(-Inf,nPar(d)),
                par.high=rep(Inf,nPar(d)),
                method='Nelder-Mead',
                control=list(fnscale=-1,maxit=10000)){
  # Bootstrap replicates
  boot.par=data.frame(matrix(NA,nsim,nPar(d)))
  names(boot.par) <- getParName(d)
  param=getParVector(d)
  if(!any(is.na(param))){
    for(i in 1:nsim){
      z=random_trim(d,n,trim,doTrunc)
      fit=estim_ML(d,z,trim,doTrunc,par.low,par.high,method,control)
      getParVector(fit)
      boot.par[i,]=getParVector(fit)
    }
  }
  return(boot.par)
}


#***************************************************************************----
# Private functions ----

#' Log-likelihood for a trimmed distribution
#'
#' @param param Numeric vector, parameters
#' @param d distributions3 object, assumed distribution.
#' @param y Numeric vector, data
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param par.low Numeric vector of size nPar(d), lower bounds for parameters.
#' @param par.high Numeric vector of size nPar(d), higher bounds for parameters.
#' @return The log-likelihood (may return NA's or -Inf)
#' @keywords internal

trimLL<-function(param,d,y,trim,doTrunc,par.low,par.high){
  if(any(is.na(param))) {return(NaN)}
  if(any(param<par.low)) {return(NaN)}
  if(any(param>par.high)) {return(NaN)}
  mask0=(y<=trim[1]);n0=sum(mask0)
  mask1=(y>=trim[2]);n1=sum(mask1)
  mask=(y>trim[1] & y<trim[2])
  n=sum(mask)
  if(n==0) {return(NaN)} # should be at least one non-trimmed value
  if(doTrunc & n!=length(y)) {return(NaN)} # all values should be within bounds for truncation
  # Make a copy of distribution object
  D=copyDist(d,param)
  # Get probabilities
  p=getPs(D,trim)
  # check feasability
  if(is.na(p$left)) {return(NaN)}
  if(p$left==1) {return(NaN)} # should be in [0,1)
  if(is.na(p$right)) {return(NaN)}
  if(p$right==1) {return(NaN)} # should be in [0,1)
  if(is.na(p$mid)) {return(NaN)}
  if(p$mid<=0) {return(NaN)} # should not be zero
  # Apply pdf to non-trimmed values
  pdfs=distributions3::log_pdf(D,y[mask])
  if(any(is.na(pdfs))) {return(NaN)}
  if(any(pdfs==-Inf)) {return(-Inf)}
  # get full likelihood
  ll=sum(pdfs)
  if(doTrunc){
    ll=ll-n*log(p$mid) # just normalize pdf
  } else {
    if(n0>0) { # at least one left-trimmed obs
      if(p$left==0) { # impossible to observe left-trimmed obs if pleft=0
        ll=-Inf
      } else {
        ll=ll+n0*log(p$left)
      }
    }
    if(n1>0) { # at least one right-trimmed obs
      if(p$right==0) { # impossible to observe right-trimmed obs if pright=0
        ll=-Inf
      } else {
        ll=ll+n1*log(p$right)
      }
    }
  }
  return(ll)
}
