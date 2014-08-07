#' 
#' @title Simulates case and controls
#' @description Generates affected and non-affected subjects until the set sample 
#' size is achieved.
#' @param n Number of observations to generate per iteration
#' @param ncases Number of cases to simulate
#' @param ncontrols Number of controls to simulate
#' @param max.sample.size Maximum number of observations allowed
#' @param pheno.prev Prevalence of the binary outcome
#' @param freq1 Minor allele frequency of the 1st genetic determinant
#' @param g1.model Genetic model of the 1st genetic determinant; 0 for binary and 1 for additive
#' @param g1.OR Odds ratios of the 1st genetic determinant
#' @param freq2 Minor allele frequency of the 2st genetic determinant
#' @param g2.model Genetic model of the 2st genetic determinant; 0 for binary and 1 for additive
#' @param g2.OR Odds ratios of the 2st genetic determinant
#' @param r Pearson coefficient of correlation for the desired level of LD
#' @param b.OR Baseline odds ratio for subject on 95 percent population 
#' centile versus 5 percentile. This parameter reflects the heterogeneity in disease 
#' risk arising from determinates that have not been measured or have not been 
#' included in the model.
#' @param ph.error misclassification rates: 1-sensitivity and 1-specificity
#' @return A matrix
#' @keywords internal
#' @author Gaye A.
#'
sim.CC.data.LD <-
function(n=NULL, ncases=NULL, ncontrols=NULL, max.sample.size=NULL, pheno.prev=NULL,
freq1=NULL, g1.model=NULL, g1.OR=NULL, freq2=NULL, g2.model=NULL, g2.OR=NULL, 
r=NULL, b.OR=NULL, ph.error=NULL)
{
   # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
   complete <- 0
   complete.absolute <- 0
   cases.complete <- 0
   controls.complete <- 0
   block <- 0

   # SET UP A MATRIX TO STORE THE GENERATED DATA
   sim.matrix <- matrix(numeric(0), ncol=7)

   # SET LOOP COUNTER
   numloops <- 0
   
   # vector to store the empirical r values
    estimated.rs <- c()
    estimated.Ds <- c()
    estimated.Dprimes <- c()
     
   # LOOP UNTIL THE SET NUMBER OF CASES AND OR CONTROLS IS ACHIEVED OR THE 
   # THE SET POPULATION SIZE TO SAMPLE FROM IS REACHED
   while(complete==0 && complete.absolute==0)
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
       estimated.rs <- append(estimated.rs, out$estimated.r)
       estimated.Ds <- append(estimated.Ds, out$estimated.D)
       estimated.Dprimes <- append(estimated.Dprimes, out$estimated.Dprime)
       LDgeno.data <- out$data
       allele.A1 <- LDgeno.data$allele.A1
       allele.B1 <- LDgeno.data$allele.B1
       allele.A2 <- LDgeno.data$allele.A2
       allele.B2 <- LDgeno.data$allele.B2
       geno1 <- LDgeno.data$geno1.U
       geno2 <- LDgeno.data$geno2.U

       # GENERATE SUBJECT EFFECT DATA THAT REFLECTS BASELINE RISK: 
       # NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR WITH APPROPRIATE 
       # VARIANCE ON SCALE OF LOG-ODDS
       s.effect.data <- sim.subject.data(n, b.OR)
					  
       # GENERATE THE TRUE OUTCOME DATA
       pheno.data <- sim.pheno.bin.LD(num.obs=n, disease.prev=pheno.prev, genotype1=geno1, genotype2=geno2, 
                                       subject.effect.data=s.effect.data, geno1.OR=g1.OR, geno2.OR=g1.OR)
       true.phenotype <- pheno.data
       
       # GENERATE THE OBSERVED OUTCOME DATA FROM WHICH WE SELECT CASES AND CONTROLS
       obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=0, pheno.error=ph.error)
       pheno <- obs.phenotype
       
       # STORE THE TRUE OUTCOME, GENETIC AND ENVIRONMENT AND ALLELE DATA IN AN OUTPUT MATRIX 
       # WHERE EACH ROW HOLDS THE RECORDS OF ONE INDIVUDAL
       sim.matrix.temp <- cbind(pheno,geno1,allele.A1,allele.B1,geno2,allele.A2,allele.B2)

       # UPDATE THE MATRIX THAT HOLDS ALL THE DATA GENERATED SO FAR, AFTER EACH LOOP
       sim.matrix <- rbind(sim.matrix, sim.matrix.temp)

       # SELECT OUT CASES
       sim.matrix.cases <- sim.matrix[pheno==1,]

       # SELECT OUT CONTROLS
       sim.matrix.controls <- sim.matrix[pheno==0,]

       # COUNT THE NUMBER OF CASES AND CONTROLS THAT HAS BEEN GENERATED
       cases.simulated <- dim(sim.matrix.cases)[1]
       controls.simulated <- dim(sim.matrix.controls)[1]

       # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
       if(cases.simulated >= ncases)
       {
         sim.matrix.cases <- sim.matrix.cases[1:ncases,]
         cases.complete <- 1
       }

       # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
       if(controls.simulated>=ncontrols)
       {
         sim.matrix.controls <- sim.matrix.controls[1:ncontrols,]
         controls.complete <- 1
       }

       # HAVE WE NOW GENERATED THE SET NUMBER OF CASES AND CONTROLS?
       complete <- cases.complete*controls.complete		

       # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
       complete.absolute <- (((block+1)*n)>=max.sample.size)
       if(complete.absolute==1) {sample.size.excess <- 1}else{sample.size.excess <- 0}

        # INCREMENT LOOP COUNTER
        numloops <- numloops + 1
       
   }

   # STACK FINAL DATA MATRIX WITH CASES FIRST
   sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # NAME THE COLUMNS OF THE MATRIX AND RETURN IT AS A DATAFRAMEDATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "genotype1", "allele.A1", "allele.B1", "genotype2", "allele.A2", "allele.B2")
   mm <- list(data=data.frame(sim.matrix), allowed.sample.size.exceeded=sample.size.excess, 
              estimated.r.value=mean(estimated.rs,na.rm=TRUE), estimated.D.value=mean(estimated.Ds,na.rm=TRUE), 
              estimated.Dprime.value=mean(estimated.Dprimes,na.rm=TRUE))
}

