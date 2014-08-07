#' 
#' @title Generates genotypes for two genetic variants in LD
#' @description  Generates alleles of two SNPs in LD and uses these alleles to form the genotypes 
#' of the genetic variants. Each variant can be binary or additive
#' @param num.obs Number of observations to simulate
#' @param MAF Minor allele frequencies of the two variants
#' @param is.add Models of the two variants
#' @param r.target Pearson coefficient of correlation for the desired level of LD
#' @param cov.mat.req Covariance matrix required required to achieved the desired LD
#' @return A dataframe that holds the genotypes and allelesdata of the two SNPs
#' @keywords internal
#' @author Gaye A.
#' @references Montana, G. 2005, \code{HapSim: a simulation tool for 
#' generating haplotype data with pre-specified allele frequencies 
#' and LD coefficients.}, Bioinformatics, \bold{vol. 21 (23)}, 
#' pp.4309-4311.
#' 
sim.LDgeno.data <-
function(num.obs=20000, MAF=c(0.1,0.1), is.add=c(0,0), r.target=0.7, cov.mat.req)
{
     maf.snp1 <- MAF[1]
     maf.snp2 <- MAF[2]
     is.add.gene <- is.add

     # CORRECTION TERM FOR MEAN CENTERING FOR ADDITIVE AND BINARY GENE
     mean.geno.add.gene <- c(2*MAF[1]*(1-MAF[1])+2*(MAF[1]^2), 2*MAF[2]*(1-MAF[2])+2*(MAF[2]^2))
     mean.geno.binary <- c(2*MAF[1]*(1-MAF[1])+(MAF[1]^2), 2*MAF[2]*(1-MAF[2])+(MAF[2]^2))

     # GENERATE ALLELES
     out1 <- sim.LDsnps (num.obs, maf.snp1, maf.snp2, r.target, cov.mat.req)
     out2 <- sim.LDsnps (num.obs, maf.snp1, maf.snp2, r.target, cov.mat.req)
     A.alleles <- out1$alleles
     B.alleles <- out2$alleles
     estimated.r <- out1$estimated.r
     estimated.D <- out1$estimated.D
     estimated.Dprime <- out1$estimated.Dprime

     # GENE 1
     allele.A1 <- A.alleles[,1]
     allele.B1 <- B.alleles[,1]
     geno1 <- allele.A1+allele.B1

     if(is.add.gene[1]==1)
     {
        geno1.U <- geno1
        geno1.U <- geno1.U-mean.geno.add.gene[1]
     }else{
        geno1.U <- geno1>0
        geno1.U <- geno1.U-mean.geno.binary[1]
     }

     # GENE 2
     allele.A2 <- A.alleles[,2]
     allele.B2 <- B.alleles[,2]     
     geno2 <- allele.A2+allele.B2
     if(is.add.gene[2]==1)
     {
        geno2.U <- geno2
        geno2.U <- geno2.U-mean.geno.add.gene[2]
     }else{
        geno2.U <- geno2 > 0
        geno2.U <- geno2.U-mean.geno.binary[2]
     }

   # RETURN A DATAFRAME
   return(list('data'=data.frame(allele.A1, allele.B1, geno1.U, allele.A2, allele.B2, geno2.U),'estimated.r'=estimated.r, 'estimated.D'=estimated.D, 'estimated.Dprime'=estimated.Dprime))
 }

