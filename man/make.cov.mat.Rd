\name{make.cov.mat}
\alias{make.cov.mat}
\title{Generates the covariance matrix required to achieved the desired LD}
\usage{
make.cov.mat(cor.mat, freqs)
}
\arguments{
  \item{cor.mat}{Correlation matrix}

  \item{freqs}{Major allele frequencies of the two snps}

  \item{is.add}{Models of the two variants}
}
\value{
A 2X2 covariance matrix
}
\description{
Finds the covariance values required to achieve the
specified frequency of major allele haplotype
}
\details{

}
\author{
Gaye A.
}
\references{
Montana, G. 2005, \code{HapSim: a simulation tool for
generating haplotype data with pre-specified allele
frequencies and LD coefficients.}, Bioinformatics,
\bold{vol. 21 (23)}, pp.4309-4311.
}
\keyword{internal}

