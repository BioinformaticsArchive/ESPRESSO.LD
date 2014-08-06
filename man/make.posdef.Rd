\name{make.posdef}
\alias{make.posdef}
\title{Turns a matrix into a positive definite one}
\usage{
make.posdef(matrix, tolerance = 1e-06)
}
\arguments{
  \item{matrix}{input matrix}

  \item{tolerance}{a constant}
}
\value{
A positive-definite matrix
}
\description{
Computes the nearest positive definite matrix of a real
symmetric matrix
}
\details{

}
\author{
Gaye A.
}
\references{
N.J. Higham, 1988 \code{Computing a nearest symmetric
positive semidefinite matrix}, Linear Algebra Appl.
\bold{vol. 103}, pp.103 118
}
\keyword{internal}

