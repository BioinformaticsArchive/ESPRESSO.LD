\name{get.critical.results.LD}
\alias{get.critical.results.LD}
\title{Summarizes the main results}
\usage{
get.critical.results.LD(scenario = NULL, pheno.model = NULL,
  geno1.model = NULL, geno2.model = NULL, sample.sizes.required = NULL,
  power = NULL, mean.beta = NULL, LDmeasures = NULL)
}
\arguments{
  \item{scenario}{Scenario number}

  \item{pheno.model}{Type of the outcome; 0 for binary and
  2 for continuous}

  \item{geno1.model}{Genetic model of the 1st SNP; 0 for
  binary and 1 for additive}

  \item{geno2.model}{Genetic model of the 1st SNP; 0 for
  binary and 1 for additive}

  \item{sample.sizes.required}{Number of cases and controls
  or number of subjects required to achieve the desired
  power}

  \item{power}{Estimated empirical power and theoretical
  (modelled) power}

  \item{mean.beta}{Mean beta value of each of the
  determinants}

  \item{LDmeasures}{the empirical Pearson r correlation
  coefficient and Lewontin's D and D'}
}
\value{
A table containing the following items: \code{gene1.model}
Model of the first genetic determinant \code{gene2.model}
Model of the second genetic determinant
\code{number.of.cases.required1} Number of cases required
to achieve the desired power for the first SNP
\code{number.of.controls.required1} Number of controls
required to achieve the desired power for the first SNP
\code{number.of.subjects.required1} Number of subjects
required to achieve the desired power for the first SNP
\code{empirical.power1} Estimated empirical powe for the
first SNP \code{modelled.power1} Power achieved under each
model with specified sample size for the first SNP
\code{estimated.OR1} Estimated odds-ratios due to shrinkage
toward the null resulting from misclassification for the
first SNP \code{estimated.effect1} Estitmated effect size
if the outocme is continuous for the first SNP
\code{number.of.cases.required2} Number of cases required
to achieve the desired power for the second SNP
\code{number.of.controls.required2} Number of controls
required to achieve the desired power for the second SNP
\code{number.of.subjects.required2} Number of subjects
required to achieve the desired power for the second SNP
\code{empirical.power2} Estimated empirical powe for the
second SNP \code{modelled.power2} Power achieved under each
model with specified sample size for the second SNP
\code{estimated.OR2} Estimated odds-ratios due to shrinkage
toward the null resulting from misclassification for the
second SNP \code{estimated.effect2} Estitmated effect size
if the outocme is continuous for the second SNP
}
\description{
Gets the number of cases and controls or subjects and the
empirical and theoretical power under each model and prints
a summary on the screen
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

