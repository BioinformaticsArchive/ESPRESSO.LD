#' 
#' @title Carries out regression analysis
#' @description Fits a conventional unconditional logistic regression model with a binary or continuous 
#' phenotype as outcome and the genetic, environmental, interaction determinants as covariates.
#' @param pheno.model Type of outcome; 0=binary and 1=continuous
#' @param observed.data A dataframe that contains covariates and outcome data
#' @return A vector containing the beta, standard-error and z-statistic of each of the covariates
#' @keywords internal
#' @author Gaye A.
#'
glm.analysis.LD <- function(pheno.model=NULL, observed.data=NULL){

  # BINARY OUTCOME
  if(pheno.model == 0){
	   # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
	   mod.glm <- glm(phenotype ~ 1+genotype1+genotype2,family=binomial,data=observed.data)
	   mod.sum <- summary(mod.glm)
  }
  
  # QUANTITATIVE OUTCOME     
  if(pheno.model == 1){
	    # FIT A GLM FOR A GAUSSIAN OUTCOME
	    mod.glm <- glm(phenotype ~ 1+genotype1+genotype2,family=gaussian,data=observed.data)
	    mod.sum <- summary(mod.glm)     
  }
  
	geno1.beta.value <- mod.sum$coefficients[2,1]
	geno1.se.value <- mod.sum$coefficients[2,2]
	geno1.z.value <- mod.sum$coefficients[2,3]
  geno2.beta.value <- mod.sum$coefficients[3,1]
  geno2.se.value <- mod.sum$coefficients[3,2]
  geno2.z.value <- mod.sum$coefficients[3,3]
	 
  # RETURN A VECTOR
  return(list(geno1.beta=geno1.beta.value, geno1.se=geno1.se.value, geno1.z=geno1.z.value,
              geno2.beta=geno2.beta.value, geno2.se=geno2.se.value, geno2.z=geno2.z.value))
}

