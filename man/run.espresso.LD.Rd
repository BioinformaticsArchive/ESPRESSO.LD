\name{run.espresso.LD}
\alias{run.espresso.LD}
\title{Runs a full ESPRESSO analysis}
\usage{
run.espresso.LD(simulation.params = NULL, pheno.params = NULL,
  geno1.params = NULL, geno2.params = NULL, scenarios2run = 1)
}
\arguments{
  \item{simulation.params}{general parameters for the
  scenario(s) to analyse}

  \item{pheno.params}{paramaters for the outcome variables}

  \item{geno1.params}{parameters for the first genetic
  determinant}

  \item{geno2.params}{parameters for the second genetic
  determinant}

  \item{scenarios2run}{the indices of the scenarios one
  wish to analyse}
}
\value{
a summary table that contains both the input parameters and
the results of the analysis
}
\description{
This function calls the functions required to run a full
ESPRESSO analysis where the model consists of an outcome
(binary or continuous) determined by two bi-allelic SNPs
that modelled as being in LD
}
\details{

}
\examples{
{

# load the table that hold the input parameters; each of the table
# hold parameters for 4 scenarios:
# scenario 1: a binary outcome determined by two binary SNPs in LD
# scenario 2: a binary outcome determined by two additive SNPs in LD
# scenario 3: a quantitative outcome determined by two binary SNPs in LD
# scenario 4: a quantitative outcome determined by two additive SNPs in LD
data(simulation.params)
data(pheno.params)
data(geno1.params)
data(geno2.params)

# run the function for the first two scenarios, two binomial models
run.espresso.LD(simulation.params, pheno.params, geno1.params, geno2.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.LD(simulation.params, pheno.params, geno1.params, geno2.params, scenarios2run=c(3,4))

}
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

