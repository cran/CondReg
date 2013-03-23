\name{condreg}
\alias{condreg}
\title{Compute the condition number with given penalty parameter}
\usage{
  condreg(data_in, kmax)
}
\arguments{
  \item{data_in}{data in one of two formats. 1) If a list,
  elements should be an orthogonal matrix Q and a vector of
  eigenvalues L from the spectral decomposition of sample
  covariance matrix of the data. 2) If a n-by-p matrix,
  each of the n rows are the p variate random
  observations.}

  \item{kmax}{scalar regularization parameter}
}
\value{
  list of condition number regularized covariance matrix s
  and its inverse invS.
}
\description{
  Compute the condition number with given penalty parameter
}
\examples{
## True covariance matrix
\dontrun{sigma <- diag(5)}
\dontrun{sigma[3,2] <- sigma[2,3] <- 0.8}

## Generate normal random samples
\dontrun{library(MASS)}
\dontrun{X <- mvrnorm(200,rep(0,5),sigma)}

## load pre-generated data
data(condreg_data)

## Covariance estimation
crcov <- condreg(X,3)

## Inspect output
str(crcov)              ## returned object
sigma.hat <- crcov$S    ## estimate of sigma matrix
omega.hat <- crcov$invS ## estimate of inverse of sigma matrix
}

