\name{sim.QTL.data.LD}
\alias{sim.QTL.data.LD}
\title{Simulates subjects for continuous outcome}
\usage{
sim.QTL.data.LD(n = NULL, ph.mean = NULL, ph.sd = NULL, freq1 = NULL,
  freq2 = NULL, g1.model = NULL, g1.efkt = NULL, g2.model = NULL,
  g2.efkt = NULL, ld = NULL, r = NULL, pheno.rel = NULL)
}
\arguments{
  \item{n}{Number of subjects to simulate}

  \item{ph.mean}{statistical mean}

  \item{ph.sd}{standard deviation}

  \item{freq1}{Minor allele frequency if the 1st of the two
  genetic variant}

  \item{freq2}{Minor allele frequency if the 2nd of the two
  genetic variant}

  \item{g1.model}{Genetic model; 0 for binary and 1 for
  additive}

  \item{g1.efkt}{Effects of the 1st genetic variant}

  \item{g2.model}{Genetic model; 0 for binary and 1 for
  additive}

  \item{g2.efkt}{Effects of the 2nd genetic variant}

  \item{r}{desired level of correlation if the two SNPs are
  modelled as being in LD.}

  \item{pheno.rel}{reliability of the assessment for a
  quantitative outcome.}
}
\value{
A matrix
}
\description{
Generates the specified number of subjects
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

