\name{sim.LDsnps}
\alias{sim.LDsnps}
\title{Generates alleles for two biallelic SNPs in Linkage Disequilibrium}
\usage{
sim.LDsnps(num.obs, maf.snp1 = 0.1, maf.snp2 = 0.1, r = 0.7, cov.mat.req)
}
\arguments{
  \item{num.obs}{Number of observations to simulate}

  \item{maf.snp1}{Minor allele frequency of the first snp}

  \item{maf.snp2}{Minor allele frequency of the second snp}

  \item{r}{Pearson coefficient of correlation for desired
  level of LD}

  \item{cov.mat.req}{Covariance matrix required required to
  achieved the desired LD}
}
\value{
A dataframe of two variables where the rows represent
haplotypes
}
\description{
Generates alleles data for pre-spcified alleles
frequencies. The covariance matrix required to achieved the
desired LD is computed and used to produce a random vector
from a bivariate normal distribution.
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

