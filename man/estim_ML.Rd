% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{estim_ML}
\alias{estim_ML}
\title{Estimate a trimmed distribution}
\usage{
estim_ML(
  d,
  y,
  trim = c(-Inf, Inf),
  doTrunc = FALSE,
  par.low = rep(-Inf, nPar(d)),
  par.high = rep(Inf, nPar(d)),
  method = "Nelder-Mead",
  control = list(fnscale = -1, maxit = 10000)
)
}
\arguments{
\item{d}{distributions3 object, assumed distribution.
The parameters specified in d are interpreted as initial guesses
and should be both feasible (e.g. no negative scale parameter)
and data-compatible (i.e. not lead to a data point having a null pdf)}

\item{y}{Numeric vector, data}

\item{trim}{Numeric vector of size 2, trimming bounds.}

\item{doTrunc}{Logical, do truncation? default FALSE, i.e. rectification is used.}

\item{par.low}{Numeric vector of size nPar(d), lower bounds for parameters.}

\item{par.high}{Numeric vector of size nPar(d), higher bounds for parameters.}

\item{method}{Character string, optimization method, see ?optim.}

\item{control}{List, controls for optimization algorithm, see ?optim.}
}
\value{
the estimated distribution as a distributions3 object.
    Note that parameters may be NA's if the estimation failed.
}
\description{
Estimate the parameters of a trimmed distribution using
maximum likelihood.
}
\examples{
# EXAMPLE 1
# Define Normal distribution
norm <- Normal(mu=1,sigma=2)
# generate data and reset all negative values to zero
y <- random(norm,200); y[y<0] <- 0
# Estimate basic / rectified / truncated Gaussian distributions.
# The latter returns NA's because y=0 values occur,
# which is not compatible with a zero-truncated distribution.
estim_ML(d=norm,y=y)
estim_ML(d=norm,y=y,trim=c(0,Inf))
estim_ML(d=norm,y=y,trim=c(0,Inf),doTrunc=TRUE)

# EXAMPLE 2
# Define Normal distribution
norm <- Normal(mu=0.75,sigma=0.25)
# generate data and remove all values outside [0;1]
z <- random(norm,200); y <- z[z>0 & z<1]
# Estimate basic / rectified / truncated Gaussian distributions.
estim_ML(d=norm,y=y)
estim_ML(d=norm,y=y,trim=c(0,1))
estim_ML(d=norm,y=y,trim=c(0,1),doTrunc=TRUE)
# Good practice: explicitly set parameter bounds, here sigma>0
estim_ML(d=norm,y=y,par.low=c(-Inf,0))
}
