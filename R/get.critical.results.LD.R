#'
#' @title Summarizes the main results
#' @description Gets the number of cases and controls or subjects and the empirical and theoretical power under each model and prints a summary on the screen
#' @param scenario Scenario number
#' @param pheno.model Type of the outcome; 0 for binary and 2 for continuous
#' @param geno1.model Genetic model of the 1st SNP; 0 for binary and 1 for additive
#' @param geno2.model Genetic model of the 1st SNP; 0 for binary and 1 for additive
#' @param sample.sizes.required Number of cases and controls or number of subjects required to achieve the desired power
#' @param power Estimated empirical power and theoretical (modelled) power
#' @param mean.beta Mean beta value of each of the determinants
#' @param LDmeasures the empirical Pearson r correlation coefficient and Lewontin's D and D'
#' @return A table containing the following items:
#' \code{gene1.model} Model of the first genetic determinant
#' \code{gene2.model} Model of the second genetic determinant
#' \code{number.of.cases.required1} Number of cases required to achieve the desired power for the first SNP
#' \code{number.of.controls.required1} Number of controls required to achieve the desired power for the first SNP
#' \code{number.of.subjects.required1} Number of subjects required to achieve the desired power for the first SNP
#' \code{empirical.power1} Estimated empirical powe for the first SNP
#' \code{modelled.power1} Power achieved under each model with specified sample size  for the first SNP
#' \code{estimated.OR1} Estimated odds-ratios due to shrinkage toward the null resulting from misclassification for the first SNP
#' \code{estimated.effect1} Estitmated effect size if the outocme is continuous for the first SNP
#' \code{number.of.cases.required2} Number of cases required to achieve the desired power for the second SNP
#' \code{number.of.controls.required2} Number of controls required to achieve the desired power for the second SNP
#' \code{number.of.subjects.required2} Number of subjects required to achieve the desired power for the second SNP
#' \code{empirical.power2} Estimated empirical powe for the second SNP
#' \code{modelled.power2} Power achieved under each model with specified sample size  for the second SNP
#' \code{estimated.OR2} Estimated odds-ratios due to shrinkage toward the null resulting from misclassification for the second SNP
#' \code{estimated.effect2} Estitmated effect size if the outocme is continuous for the second SNP
#' @keywords internal
#' @author Gaye A.
#'
get.critical.results.LD <-
function(scenario=NULL, pheno.model=NULL,geno1.model=NULL,geno2.model=NULL,
         sample.sizes.required=NULL,power=NULL,mean.beta=NULL, LDmeasures=NULL)
{
		 
 	 if(geno1.model==0){
 	   g1.model <- "binary"
 	 }else{
 	   g1.model <- "additive"
 	 }
    
 	 if(geno2.model==0){
 	   g2.model <- "binary"
 	 }else{
 	   g2.model <- "additive"
 	 }

   if(pheno.model == 0){
			numcases1 <- sample.sizes.required[[1]]
			numcontrols1 <- sample.sizes.required[[2]]
			numcases2 <- sample.sizes.required[[3]]
			numcontrols2 <- sample.sizes.required[[4]]
   }else{
			numsubjects1 <- sample.sizes.required[[1]]
			numsubjects2 <- sample.sizes.required[[2]]
   }

  if(pheno.model==0){
     # estimated ORs
     estimated.OR1 <- exp(mean.beta[1])
     estimated.effect1 <- 'NA'
     estimated.OR2 <- exp(mean.beta[2])
     estimated.effect2 <- 'NA'

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat(" Outcome: binary \n")
					cat(" First genetic determinant:",g1.model,"\n")
					cat(" Second genetic determinant:",g2.model)

					cat("\n\nNumber of cases required\n")
					cat("------------------------\n")
					cat(" First genetic determinant:", numcases1,"\n")
          cat(" Second genetic determinant:", numcases2)

					cat("\n\nNumber of controls required\n")
					cat("---------------------------\n")
          cat(" First genetic determinant:", numcontrols1,"\n")
          cat(" Second genetic determinant:", numcontrols2)
     
  				cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" First genetic determinant:",power[[1]], "\n")
          cat(" Second genetic determinant:",power[[3]])

					cat("\n\nModel power\n")
					cat("-----------\n")
					cat(" First genetic determinant:",round(power[[2]],2), "\n")
          cat(" Second genetic determinant:",round(power[[4]],2))

					cat("\n\nEstimated effect size\n")
					cat("-----------\n")
					cat(" First genetic determinant:",round(estimated.OR1,2), "\n")
          cat(" Second genetic determinant:",round(estimated.OR2,2))
     
          cat("\n\nEstimated LD\n")
          cat("-----------\n")
          cat(" Pearson r:",LDmeasures[1], "\n")
          cat(" Lewontin D:",LDmeasures[2], "\n")
          cat(" Lewontin D':",LDmeasures[3])

					cat("\n\n---- END OF SUMMARY ----\n")

		   crit.res <- c(g1.model,g2.model,numcases1,numcontrols1,round(power[[1]],2),round(power[[2]],2),
		                 round(estimated.OR1,2),estimated.effect1,numcases2,numcontrols2,round(power[[2]],2),
                     round(power[[2]],2),round(estimated.OR2,2), estimated.effect2)
     
		   return(list(gene1.model=crit.res[1], gene2.model=crit.res[2],number.of.cases.required1=crit.res[3],
		               number.of.controls.required1=crit.res[4],empirical.power1=crit.res[5], 
		               modelled.power1=crit.res[6],estimated.OR1=crit.res[7], estimated.effect1=crit.res[8],
		               number.of.cases.required2=crit.res[9], number.of.controls.required2=crit.res[10],empirical.power2=crit.res[11], 
		               modelled.power2=crit.res[12],estimated.OR2=crit.res[13], estimated.effect2=crit.res[14]))

  }else{
     # estimated ORs
     estimated.effect1 <- mean.beta[1]
     estimated.OR1 <- 'NA'
     estimated.effect2 <- mean.beta[2]
     estimated.OR2 <- 'NA'

         cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
         cat("\nModels\n")
         cat("------\n")
         cat(" Outcome: binary \n")
         cat(" First genetic determinant:",g1.model,"\n")
         cat(" Second genetic determinant:",g2.model)
         
         cat("\n\nNumber of subjects required\n")
         cat("------------------------\n")
         cat(" First genetic determinant:", numsubjects1,"\n")
         cat(" Second genetic determinant:", numsubjects2)

         cat("\n\nEmpirical power\n")
         cat("---------------\n")
         cat(" First genetic determinant:",power[[1]], "\n")
         cat(" Second genetic determinant:",power[[3]])
         
         cat("\n\nModel power\n")
         cat("-----------\n")
         cat(" First genetic determinant:",round(power[[2]],2), "\n")
         cat(" Second genetic determinant:",round(power[[4]],2))
         
         
         cat("\n\nEstimated effect size\n")
         cat("-----------\n")
         cat(" First genetic determinant:",round(estimated.effect1,2), "\n")
         cat(" Second genetic determinant:",round(estimated.effect2,2))
     
         cat("\n\nEstimated LD\n")
         cat("-----------\n")
         cat(" Pearson r:",LDmeasures[1], "\n")
         cat(" Lewontin D:",LDmeasures[2], "\n")
         cat(" Lewontin D':",LDmeasures[3])

				 cat("\n\n---- END OF SUMMARY ----\n")

		   crit.res <- c(g1.model,g2.model,numsubjects1,round(power[[1]],2),round(power[[2]],2),estimated.OR1,round(estimated.effect1,2),
		                 numsubjects2,round(power[[3]],2),round(power[[4]],2),estimated.OR2,round(estimated.effect2,2))
		   return(list(gene1.model=crit.res[1], gene2.model=crit.res[2],number.of.subjects.required1=crit.res[3],
		              empirical.power1=crit.res[4], modelled.power1=crit.res[5], estimated.OR1=crit.res[6],estimated.effect1=crit.res[7],
		              number.of.subjects.required2=crit.res[8],empirical.power2=crit.res[9], modelled.power2=crit.res[10], 
                  estimated.OR2=crit.res[11],estimated.effect2=crit.res[12]))
  }
}

