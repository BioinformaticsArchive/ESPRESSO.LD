\name{sim.LDgeno.data}
\alias{sim.LDgeno.data}
\title{Generates genotypes for two genetic variants in LD}
\usage{
sim.LDgeno.data(num.obs = 20000, MAF = c(0.1, 0.1), is.add = c(0, 0),
  r.target = 0.7, cov.mat.req)
}
\arguments{
  \item{num.obs}{Number of observations to simulate}

  \item{MAF}{Minor allele frequencies of the two variants}

  \item{is.add}{Models of the two variants}

  \item{r.target}{Pearson coefficient of correlation for
  desired level of LD}

  \item{cov.mat.req}{Covariance matrix required required to
  achieved the desired LD}
}
\value{
A dataframe that holds the genotypes and allelesdata of the
two SNPs
}
\description{
Generates alleles of two SNPs in LD and uses these alleles
to form the genotypes of the genetic variants. Each variant
can be binary or additive
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

