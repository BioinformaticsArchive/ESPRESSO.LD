#' 
#' @title Calculates the empirical and theoretical power
#' @description The function determines the empirical and theoretical power. 
#' The empirical power is the proportion of simulations in which 
#' the z-statistic for the parameter of interest exceeds the z-statistic 
#' for the desured level if statistical significance. 
#' The theoretical power is the power of the study.
#' @param pval cut-off p-value defining statistical significance.
#' @param z.values z-statistic of the determinant.
#' @param mean.model.z mean z-statistic values of the two SNPs
#' @return a list that contains the computed empirical power and theoretical power.
#' @keywords internal
#' @author Gaye A.
#'
power.calc.LD <- function(pval=NULL, z.values=NULL, mean.model.z=NULL){
  
  if(is.null(z.values)){
    message("ALERT!")
    message(" No z-statistics found")
    cat(" Check the argument 'z.values'")
    stop("End of process!", call.=FALSE)
  }
  
  if(is.null(mean.model.z)){
    message("ALERT!")
    message(" The argument 'mean.model.z' is set to NULL.")
    message(" This argument should be the ratio 'mean.beta/mean.se'.")
    stop("End of process!", call.=FALSE)
  }
  
  # CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE 
  z.pval <- qnorm(1-pval/2)
  
  # GET EMPIRICAL POWER: THE PROPORTION OF SIMULATIONS IN WHICH THE 
  # Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC 
  # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
  geno1.empirical.power <- round(mean((z.values[[1]] > z.pval), na.rm=TRUE),3)
  geno2.empirical.power <- round(mean((z.values[[2]] > z.pval), na.rm=TRUE),3)
  
  # GET THE MODELLED POWER
  geno1.modelled.power <- pnorm(mean.model.z[1]-z.pval)
  geno2.modelled.power <- pnorm(mean.model.z[2]-z.pval)
  
  return(list(geno1.empirical=geno1.empirical.power, geno1.modelled=geno1.modelled.power,
              geno2.empirical=geno2.empirical.power, geno2.modelled=geno2.modelled.power))
}
