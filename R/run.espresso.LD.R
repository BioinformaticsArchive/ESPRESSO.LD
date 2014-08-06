#' 
#' @title Runs a full ESPRESSO analysis
#' @description This function calls the functions required to run a full ESPRESSO analysis 
#'  where the model consists of an outcome (binary or continuous) determined by two bi-allelic SNPs that
#' modelled as being in LD
#' @param simulation.params general parameters for the scenario(s) to analyse
#' @param pheno.params paramaters for the outcome variables
#' @param geno1.params parameters for the first genetic determinant
#' @param geno2.params parameters for the second genetic determinant
#' @param scenarios2run the indices of the scenarios one wish to analyse
#' @return a summary table that contains both the input parameters and 
#' the results of the analysis
#' @export
#' @author Gaye A.
#' @references Montana, G. 2005, \code{HapSim: a simulation tool for 
#' generating haplotype data with pre-specified allele frequencies 
#' and LD coefficients.}, Bioinformatics, \bold{vol. 21 (23)}, 
#' pp.4309-4311.
#' @examples {
#'   
#' # load the table that hold the input parameters; each of the table
#' # hold parameters for 4 scenarios:
#' # scenario 1: a binary outcome determined by two binary SNPs in LD
#' # scenario 2: a binary outcome determined by two additive SNPs in LD
#' # scenario 3: a quantitative outcome determined by two binary SNPs not in LD
#' # scenario 4: a quantitative outcome determined by two additive SNPs not in LD 
#' data(simulation.params) 
#' data(pheno.params)
#' data(geno1.params)
#' data(geno2.params)
#' 
#' # run the function for the first two scenarios, two binomial models
#' run.espresso.LD(simulation.params, pheno.params, geno1.params, geno2.params, scenarios2run=c(1,2))
#'
#' # run the function for the last two scenarios, two gaussian models
#' run.espresso.LD(simulation.params, pheno.params, geno1.params, geno2.params, scenarios2run=c(3,4))
#' 
#' }
#'
run.espresso.LD <- function(simulation.params=NULL, pheno.params=NULL, geno1.params=NULL, 
                             geno2.params=NULL, scenarios2run=1){

# IF AN INPUT FILE IS NOT SUPPLIED LOAD THE DEFAULT TABLES WARNING
if(is.null(simulation.params)){
  cat("\n WARNING!\n")
  cat(" No simulation parameters supplied\n")
  cat(" The default simulation parameters will be used\n")
  simulation.params <- data(simulation.params)
}

if(is.null(pheno.params)){
  cat("\n WARNING!\n")
  cat(" No outcome parameters supplied\n")
  cat(" The default outcome parameters will be used\n")
  pheno.params <- data(pheno.params)
}

if(is.null(geno1.params)){
  cat("\n WARNING!\n")
  cat(" No genotype 1 parameters supplied\n")
  cat(" The default genotype 1 parameters will be used\n")
  geno1.params <- data(geno1.params)
}

if(is.null(geno2.params)){
  cat("\n WARNING!\n")
  cat(" No genotype 2 parameters supplied\n")
  cat(" The default genotype 2 parameters will be used\n")
  geno2.params <- data(geno2.params)
}

# MERGE INPUT FILES TO MAKE ONE TABLE OF PARAMETERS
s.temp1 <- merge(simulation.params, pheno.params)
s.temp2 <- merge(s.temp1, geno1.params)
s.parameters <- merge(s.temp2, geno2.params)


#----------LOAD SET UP UP INITIAL PARAMETERS------------#

# PRINT TRACER CODE EVERY Nth ITERATION
# THIS ENSURES THAT YOU CAN SEE IF THE PROGRAM GRINDS TO A HALT FOR SOME REASON (IT SHOULDN'T)
trace.interval <- 1


# CREATE UP TO 20M SUBJECTS IN BLOCKS OF 20K UNTIL REQUIRED NUMBER OF
# CASES AND CONTROLS IS ACHIEVED. IN GENERAL THE ONLY PROBLEM IN ACHIEVING THE
# REQUIRED NUMBER OF CASES WILL OCCUR IF THE DISEASE PREVALENCE IS VERY LOW
allowed.sample.size <- 20000000
block.size <- 20000


# DECLARE MATRIX THAT STORE THE RESULTS FOR EACH SCENARIO (ONE PER SCENARIO PER ROW)
output.file <- "output.csv"
output.matrix <- matrix(numeric(0), ncol=48)
column.names <- c(colnames(s.parameters), "exceeded.sample.size?", "numcases.required.geno1", "numcontrols.required.geno1", 
                 "numsubjects.required.geno1", "empirical.power.geno1", "modelled.power.geno1", "estimated.OR.geno1", 
                 "estimated.effect.geno1", "numcases.required.geno2", "numcontrols.required.geno2", 
                 "numsubjects.required.geno2", "empirical.power.geno2", "modelled.power.geno2", "estimated.OR.geno2", 
                 "estimated.effect.geno2", "estimated.r", "estimated.Lewontin.D", "estimated.Lewontin.D'")
write(t(column.names),output.file,dim(output.matrix)[2],append=TRUE,sep=";")
                 
#-----------LOOP THROUGH THE SCENARIOS - DEALS WITH ONE SCENARIO AT A TIME-------------

for(j in c(scenarios2run))
{

   # RANDOM NUMBER GENERATOR STARTS WITH SEED SET AS SPECIFIED 
   set.seed(s.parameters$seed.val[j])

   # SIMULATION PARAMETERS
   scenario.id <- s.parameters$scenario.id[j]         
   seed.val <- s.parameters$seed.val[j]               
   numsims <- s.parameters$numsims[j]                 
   numcases <- s.parameters$numcases[j]               
   numcontrols <- s.parameters$numcontrols[j]  
   numsubjects <- s.parameters$numsubjects[j]
   LD <- s.parameters$LD[j]
   r.target <- s.parameters$r.target[j]            
   baseline.OR <- s.parameters$RR.5.95[j]                
   pval <- s.parameters$p.val[j]                      
   power <- s.parameters$power[j]
   
   # OUTCOME PARAMETERS
   pheno.model <- s.parameters$pheno.model[j]
   pheno.mean <- s.parameters$pheno.mean[j]
   pheno.sd <- s.parameters$pheno.sd[j]
   disease.prev <- s.parameters$disease.prev[j]
   pheno.error <- c(1-s.parameters$pheno.sensitivity[j],1-s.parameters$pheno.specificity[j])
   pheno.reliability <- s.parameters$pheno.reliability[j]    

   # GENETIC DETERMINANT 1 PARAMETERS
   geno1.model<- s.parameters$geno1.model[j]
   MAF1 <-  s.parameters$MAF1[j]    
   geno1.OR <- s.parameters$geno1.OR[j]
   geno1.efkt <- s.parameters$geno1.efkt[j]
   geno1.error <- c(1-s.parameters$geno1.sensitivity[j],1-s.parameters$geno1.specificity[j])
   
   # GENETIC DETERMINANT 2 PARAMETERS
   geno2.model<- s.parameters$geno2.model[j]
   MAF2 <-  s.parameters$MAF2[j]    
   geno2.OR <- s.parameters$geno2.OR[j]
   geno2.efkt <- s.parameters$geno2.efkt[j]
   geno2.error <- c(1-s.parameters$geno2.sensitivity[j],1-s.parameters$geno2.specificity[j])
   

   # VECTORS TO HOLD BETA, SE AND Z VALUES AFTER EACH RUN OF THE SIMULATION
   geno1.beta.values <- rep(NA,numsims)
   geno1.se.values <- rep(NA,numsims)
   geno1.z.values <-rep(NA,numsims)
   geno2.beta.values <- rep(NA,numsims)
   geno2.se.values <- rep(NA,numsims)
   geno2.z.values <- rep(NA,numsims)
   estimated.rs <- rep(NA,numsims)
   estimated.Ds <- rep(NA,numsims)
   estimated.Dprimes <- rep(NA,numsims)

   # TRACER TO DETECT EXCEEDING MAX ALLOWABLE SAMPLE SIZE
   sample.size.excess <- 0

   # GENERATE AND ANALYSE DATASETS ONE AT A TIME 
   for(s in 1:numsims)            # s from 1 to total number of simulations
   {

      #----------------------------------GENERATE "TRUE" DATA-----------------------------#
      

      if(pheno.model == 0){ # UNDER BINARY OUTCOME MODEL
        # GENERATE CASES AND CONTROLS UNTILL THE REQUIRED NUMBER OF CASES, CONTROLS IS ACHIEVED 
        sim.data <- sim.CC.data.LD(n=block.size, ncases=numcases, ncontrols=numcontrols, 
                                    max.sample.size=allowed.sample.size, pheno.prev=disease.prev,
                                    freq1=MAF1, g1.model=geno1.model, g1.OR=geno1.OR, 
                                    freq2=MAF2, g2.model=geno2.model, g2.OR=geno2.OR, 
                                    ld=LD, r=r.target, b.OR=baseline.OR, ph.error=pheno.error)
      }else{ # UNDER QUANTITATIVE OUTCOME MODEL
        # GENERATE THE SPECIFIED NUMBER OF SUBJECTS
        sim.data <- sim.QTL.data.LD(n=numsubjects,ph.mean=pheno.mean,ph.sd=pheno.sd,freq1=MAF1,
                                      g1.model=geno1.model,g1.efkt=geno1.efkt,freq2=MAF2,
                                      g2.model=geno2.model,g2.efkt=geno2.efkt,ld=LD, r=r.target,
                                      pheno.rel=pheno.reliability)
      }
      true.data <- sim.data$data
      
      if(LD == 1){
        estimated.rs[s] <- sim.data$estimated.r.value
        estimated.Ds[s] <- sim.data$estimated.r.value
        estimated.Dprimes[s] <- sim.data$estimated.r.value
      }       

      #------------SIMULATE ERRORS AND ADD THEM TO THE TRUE COVARIATES DATA TO OBTAIN OBSERVED COVARIATES DATA-----------#

      # ADD APPROPRIATE ERRORS TO PRODUCE OBSERVED GENOTYPES 
      observed.data <- get.observed.data.LD(data=true.data,g1.error=geno1.error,g1.model=geno1.model,freq1=MAF1,
                                             g2.error=geno2.error,g2.model=geno2.model,freq2=MAF2)

      #--------------------------DATA ANALYSIS ----------------------------#

      glm.estimates <- glm.analysis.LD(pheno.model, observed.data)
     
      geno1.beta.values[s] <- glm.estimates$geno1.beta
      geno1.se.values[s] <- glm.estimates$geno1.se
      geno1.z.values[s] <- glm.estimates$geno1.z
      geno2.beta.values[s] <- glm.estimates$geno2.beta
      geno2.se.values[s] <- glm.estimates$geno2.se
      geno2.z.values[s] <- glm.estimates$geno2.z
      
      # PRINT TRACER AFTER EVERY Nth DATASET CREATED
      message(s," of ",numsims," runs completed in scenario ",scenario.id)

   }
   cat("\n\n")

   #------------------------ SUMMARISE RESULTS ACROSS ALL SIMULATIONS---------------------------#

   # SUMMARISE PRIMARY PARAMETER ESTIMATES
   # COEFFICIENTS ON LOG-ODDS SCALE
   geno1.mean.beta <- round(mean(geno1.beta.values, na.rm=T),3)
   geno1.mean.se <- round(sqrt(mean(geno1.se.values^2, na.rm=T)),3)
   geno1.mean.model.z <- geno1.mean.beta/geno1.mean.se
   
   geno2.mean.beta <- round(mean(geno2.beta.values, na.rm=T),3)
   geno2.mean.se <- round(sqrt(mean(geno2.se.values^2, na.rm=T)),3)
   geno2.mean.model.z <- geno2.mean.beta/geno2.mean.se
   
 
   #---------------------------POWER AND SAMPLE SIZE CALCULATIONS----------------------#
   mean.model.z <- c(geno1.z.values, geno2.z.values)
   z.values <- list(geno1.mean.model.z, geno2.mean.model.z)
   
   # CALCULATE THE SAMPLE SIZE REQUIRED UNDER EACH MODEL
   sample.sizes.required <- samplsize.calc.LD(numcases, numcontrols, numsubjects, pheno.model, pval, power, mean.model.z)

   # CALCULATE EMPIRICAL POWER AND THE MODELLED POWER 
   # THE EMPIRICAL POWER IS SIMPLY THE PROPORTION OF SIMULATIONS IN WHICH
   # THE Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC
   # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
   power <- power.calc.LD(pval, z.values, mean.model.z)


   #------------------MAKE FINAL A TABLE THAT HOLDS BOTH INPUT PARAMETERS AND OUTPUT RESULTS---------------#
   mean.beta <- c(geno1.mean.beta, geno2.mean.beta)
   if(LD == 0){estimated.r <- NA}else{estimated.r <- round(mean(estimated.rs, na.rm=TRUE),2)}
   if(LD == 0){estimated.D <- NA}else{estimated.D <- round(mean(estimated.Ds, na.rm=TRUE),2)}
   if(LD == 0){estimated.Dprime <- NA}else{estimated.Dprime <- round(mean(estimated.Dprimes, na.rm=TRUE),2)}
   estimatedLDcoeffs <- c(estimated.r, estimated.D, estimated.Dprime)
   
   critical.res <- get.critical.results.LD(j,pheno.model,geno1.model,geno2.model,sample.sizes.required,power,mean.beta,estimatedLDcoeffs)

   #  WHEN OUTCOME IS BINARY INFORM IF RECORD EXCEEDED MAXIMUM SAMPLE SIZE
   if(pheno.model==0){
     sample.size.excess <- sim.data$allowed.sample.size.exceeded
     if(sample.size.excess==1)
     {
       excess <- "yes"
       cat("\nTO GENERATE THE NUMBER OF CASES SPECIFIED AT OUTSET\n")
       cat("THE SIMULATION EXCEEDED THE MAXIMUM POPULATION SIZE OF ", allowed.sample.size,"\n")
     }else{
       excess <- "no"
     }
   }
   
   inparams <- s.parameters[j,]
   if(LD==0){ s.parameters[,8] <- 'NA'}
   if(pheno.model==0){
      mod <- "binary"
      inparams [c(6,14,15,18,22,28)] <- "NA"
      inputs <- inparams
      outputs <- c(excess, critical.res[[3]], critical.res[[4]], "NA", critical.res[[5]], 
                   critical.res[[6]], critical.res[[7]], critical.res[[8]], critical.res[[9]],
                   critical.res[[10]], "NA", critical.res[[11]], critical.res[[12]], 
                   critical.res[[13]], critical.res[[14]])
   }else{
      mod <- "quantitative"
      inparams [c(4,5,9,13,23,24,29,30)] <- "NA"
      inputs <- inparams
      outputs <- c("NA", "NA", "NA", critical.res[[3]], critical.res[[4]], critical.res[[5]], 
                   critical.res[[6]], critical.res[[7]], "NA", "NA", critical.res[[8]], 
                   critical.res[[9]], critical.res[[10]],critical.res[[11]], critical.res[[12]])
   }
   
   jth.row <- as.character(c(inputs,outputs,estimated.r,estimated.D,estimated.Dprime))
   write(t(jth.row),output.file,dim(output.matrix)[2],append=TRUE,sep=";")
}
}


