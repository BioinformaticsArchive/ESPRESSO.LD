#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param n Number of subjects to simulate
#' @param ph.mean statistical mean
#' @param ph.sd standard deviation
#' @param freq1 Minor allele frequency if the 1st of the two genetic variant
#' @param freq2 Minor allele frequency if the 2nd of the two genetic variant
#' @param g1.model Genetic model; 0 for binary and 1 for additive
#' @param g1.efkt Effects of the 1st genetic variant
#' @param g2.model Genetic model; 0 for binary and 1 for additive
#' @param g2.efkt Effects of the 2nd genetic variant
#' @param r Pearson coefficient of correlation for the desired level of LD.
#' @param pheno.rel reliability of the assessment for a quantitative outcome.
#' @return A matrix
#' @keywords internal
#' @author Gaye A.
#'
sim.QTL.data.LD <-
function(n=NULL,ph.mean=NULL,ph.sd=NULL,freq1=NULL,freq2=NULL,g1.model=NULL,g1.efkt=NULL,
         g2.model=NULL,g2.efkt=NULL,ld=NULL, r=NULL,pheno.rel=NULL)
{
  
   # the covariance matrix required to generate 2 variants with 
   # the desired ld
   cor.mat <- matrix(c(1,r,r,1),2,2) # cor. matrix
   cov.mat.req <- make.cov.mat(cor.mat, c(1-freq1, 1-freq2))
     
   # if the required covariance matrix is not positive-definite get 
   # the nearest positive-definite matrix (tolerance = 1e-06)
   if(!is.posdef(cov.mat.req, 0.000001)){
     cov.mat.req <- make.posdef(cov.mat.req, 0.000001)
   }
   # GENERATE THE TRUE GENOTYPE DATA FOR THE 1st and 2nd DETERMINANT 
   out <- sim.LDgeno.data(n, c(freq1,freq2), c(g1.model,g2.model), r, cov.mat.req)
   estimated.r <- out$estimated.r
   estimated.D <- out$estimated.D
   estimated.Dprime <- out$estimated.Dprime
   LDgeno.data <- out$data
   allele.A1 <- LDgeno.data$allele.A1
   allele.B1 <- LDgeno.data$allele.B1
   allele.A2 <- LDgeno.data$allele.A2
   allele.B2 <- LDgeno.data$allele.B2
   geno1 <- LDgeno.data$geno1.U
   geno2 <- LDgeno.data$geno2.U
           
   # GENERATE THE TRUE OUTCOME DATA
   pheno.data <- sim.pheno.qtl.LD(numsubjects=n,pheno.mean=ph.mean,pheno.sd=ph.sd,genotype1=geno1,
                                  genotype2=geno2,geno1.efkt=g1.efkt,geno2.efkt=g2.efkt)
   true.phenotype <- pheno.data
   
   # GENERATE THE OBSERVED OUTCOME DATA 
   obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, 
                                  pheno.sd=ph.sd, pheno.reliability=pheno.rel)
   pheno <- obs.phenotype
   

   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(pheno,geno1,allele.A1,allele.B1,geno2,allele.A2,allele.B2)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id","phenotype","genotype1","allele.A1","allele.B1","genotype2","allele.A2","allele.B2")
   mm <- list(data=data.frame(sim.matrix), estimated.r.value=estimated.r, estimated.D.value=estimated.D, estimated.Dprime.value=estimated.Dprime)
}

