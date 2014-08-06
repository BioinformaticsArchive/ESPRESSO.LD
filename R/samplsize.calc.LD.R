#'  
#' @title Calculates the sample size required to achieve the desired power
#' @description Estimates by how much the simulated study size needs to be 
#' inflated or shrank in order to obtain the specified level of power. 
#' The ratio of z-statistic required for desired power to mean model z-statistic 
#' obtained indicates the relative changes in standard error required. 
#' This corresponds to the relative change on scale and square root of sample size.
#' @param numcases number of cases when outcome is binary.
#' @param numcontrols number of controls when outcome is binary.
#' @param num.subjects number of subjects when outcome is quantitative.
#' @param pheno.model outcome model, binary=0 and quantitafive=1.
#' @param pval cut-off p-value defining statistical significance.
#' @param power desired power
#' @param mean.model.z ratio of mean beta estimate over mean se estimate for each of the two SNPs
#' @return A table containing:
#' \code{numcases.required} number of cases required to achieve the desired power under binary outcome model
#' \code{numcontrols.required} number of controls required to achieve the desired power under binary outcome model
#' \code{numsubjects.required} number of subjects required to achieve the desired power under a quantitative outcome model
#' @keywords internal
#' @author Gaye A.
#'
samplsize.calc.LD <-
  function(numcases=NULL,numcontrols=NULL,num.subjects=NULL,pheno.model=NULL,pval=NULL,power=NULL,mean.model.z=NULL){
    
    # CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE AND POWER
    z.pval <- qnorm(1-pval/2)
    z.power.required <- qnorm(power)+z.pval
    
    # ESTIMATE HOW MUCH THE SIMULATED STUDY SIZE NEEDS TO BE INFLATED OR SHRINKED 
    # IN ORDER TO OBTAIN A POWER OF 80%. THE RATIO OF Z STATISTIC REQUIRED FOR DESIRED 
    # POWER TO MEAN MODEL Z STATISTIC OBTAINED INDICATES RELATIVE CHANGE REQUIRED IN STANDARD ERROR. 
    # THIS CORRESPONDS TO RELATIVE CHANGE ON SCALE OF SQUARE ROOT OF SAMPLE SIZE. RATIO OF SAMPLE 
    # SIZE IS THEREFORE THIS RATIO SQUARED.
    sample.size.inflation.required.geno1 <- (z.power.required/mean.model.z[1])^2
    sample.size.inflation.required.geno2 <- (z.power.required/mean.model.z[2])^2
    
    if(pheno.model==0){ # IF THE OUTCOME IS BINARY
      # MULTIPLY THE INPUT NUMBER OF CASES AND CONTROLS BY THIS
      # SQUARED RATIO TO GET THE REQUIRED NUMBER OF CASES AND CONTROLS
      # FOR THE DESIRED POWER
      geno1.cases <- round(numcases*sample.size.inflation.required.geno1,0)
      geno1.controls <- round(numcontrols*sample.size.inflation.required.geno1,0) 
      geno2.cases <- round(numcases*sample.size.inflation.required.geno2,0)
      geno2.controls <- round(numcontrols*sample.size.inflation.required.geno2,0) 
      return(list(geno1.numcases.required=geno1.cases, geno1.numcontrols.required=geno1.controls,
                  geno2.numcases.required=geno2.cases, geno2.numcontrols.required=geno2.controls))
      
    }else{ # IF THE OUTCOME IS CONTINUOUS
      # MULTIPLY THE INPUT NUMBER OF SUBJECTS BY THE INFLATION REQUIRED
      geno1.subjects <- round(num.subjects*sample.size.inflation.required.geno1,0) 
      geno2.subjects <- round(num.subjects*sample.size.inflation.required.geno2,0)   
      return(list(geno1.numsubjects.required=geno1.subjects, geno2.numsubjects.required=geno2.subjects))
    }
  }
