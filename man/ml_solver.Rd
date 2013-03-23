\name{ml_solver}
\alias{ml_solver}
\title{Compute shrinkage of eigenvalues for condreg}
\usage{
  ml_solver(L, k)
}
\arguments{
  \item{L}{vector of eigenvalues}

  \item{k}{vector of penalties}
}
\value{
  list of vector of shrinked eigenvalues \code{Lbar},
  optimal u value \code{uopt} and interval indicator
  \code{intv}.
}
\description{
  Compute shrinkage of eigenvalues for condreg
}

