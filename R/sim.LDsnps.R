#' 
#' @title Generates alleles for two biallelic SNPs in Linkage Disequilibrium
#' @description    Generates alleles data for pre-spcified alleles frequencies. 
#' The covariance matrix required to achieved the desired LD is computed and used to 
#' produce a random vector from a bivariate normal distribution.
#' @param num.obs Number of observations to simulate
#' @param maf.snp1 Minor allele frequency of the first snp
#' @param maf.snp2 Minor allele frequency of the second snp
#' @param r Pearson coefficient of correlation for desired level of LD
#' @param cov.mat.req Covariance matrix required required to achieved the desired LD
#' @return A dataframe of two variables where the rows represent haplotypes 
#' @keywords internal
#' @author Gaye A.
#' @references Montana, G. 2005, \code{HapSim: a simulation tool for 
#' generating haplotype data with pre-specified allele frequencies 
#' and LD coefficients.}, Bioinformatics, \bold{vol. 21 (23)}, 
#' pp.4309-4311.
#'
sim.LDsnps <-
function(num.obs, maf.snp1=0.1, maf.snp2=0.1, r=0.7, cov.mat.req)
{
			# IF ONE OF THE SNPS HAS A MAF OF 0 OR 1, DISPLAY AN ALERT AND STOP THE PROCESS
			if(maf.snp1==0 || maf.snp1==1 || maf.snp2==0 || maf.snp2==1){
			  message("ALERT!")
			  message(" Minor Allele Frequency was set to 1 or 0.")
			  message(" Lewontin's D cannnot be computed.")
			stop("End of process!", call.=FALSE)
			}

			# NUMBER OF DIALLELIC SNPS
			num.snps <- 2

			# NUMBER OF HAPLOTYPES TO GENERATE
			num.haps <- num.obs

			# TARGET MAJOR ALLELE FREQUENCIES OF THE SNPS 
			majs <- c(1-maf.snp1, 1-maf.snp2)
			freqs <- append(majs, (1-majs))

			# TARGET CORRELATION BETWEEN SNPS (LD COEFFICIENT)
			r.target <- r

			# QUANTILE VALUES CORRESPONDING TO THE SET FREQUENCIES
			quantiles <- abs(qnorm(majs))

			# TARGET LEWONTIN D
			D.target <- r.target * sqrt(prod(freqs))

			# TARGET LEWONTIN D'
			if(D.target < 0){
			   Dprime.target <- D.target / min(majs[1]*(1-majs[2]), majs[2]*(1-majs[1]))
			}else{
			   Dprime.target <- D.target / max(-majs[1]*majs[2], -(1-majs[1])*(1-majs[2]))
			}

			# FREQUENCY OF MAJOR ALLELE HAPLOTYPE (P.AB, A MAJOR ALLELE OF SNP1 AND B MAJOR ALLELE OF SNP2)
			# D = P.AB - (PA*PB) SO USING THE EQUATION D = R * SQRT(PROD(FREQS)) WE CAN COMPUTE P.AB AS BELOW
			P.AB.target <- (r.target * sqrt(prod(freqs))) + prod(majs)

			# COVARIANCE BETWEEN SNPS, WE HAVE A BERNOULI DISTRIBUTION
			# SO GET THE VARIANCES USING THE MAJS AND CALCULATE 
			# THE COVARIANCES USING THE R.TARGETS AND THE VARIANCES
			# IN BERNOULI DISRIBUTION VAR=P(1-P)
			vars <- majs*(1-majs) 
			cov.coeff <- r.target * (prod(sqrt(vars)))

			# A MATRIX TO HOLD THE GENERATED HAPLOTYPES (COLUMNS ARE SNPS AND ROWS ARE HAPLOTYPES)
			hap.matrix <- matrix(0, nrow = num.haps, ncol = num.snps)
			columns <- c()
			for(i in 1:num.snps){ columns <-  append(columns, paste("SNP", i, sep=""))}
			colnames(hap.matrix) <- columns
			rows <- c()
			for(i in 1:num.haps){ rows <-  append(rows, paste("HAP", i, sep=""))}
			row.names(hap.matrix) <- rows

			# GENERATE A MATRIX OF SAMPLES WITH A SET MULTIVARIATE NORMAL DISTRIBUTION
			temp.matrix <- mvrnorm(num.haps, mu = rep(0, num.snps), Sigma = cov.mat.req)

			# WHENEVER A VALUE OF TEMP.MATRIX >= TO THE INDEXWISE CORRESPONDING
			# QUANTILE VALUE OF A SNP, SET THE CORRESPONDING ALLELE TO 1.
			indx.snp1 <- which(temp.matrix[,1] >= quantiles[1])
			hap.matrix[indx.snp1,1] <- 1
			indx.snp2 <- which(temp.matrix[,2] >= quantiles[2])
			hap.matrix[indx.snp2,2] <- 1

			# ---- IN THE BELOW LINE 'A1' AND 'B1' ARE RESPECTIVELY THE MAJOR AND MINOR ALLELE OF SNP1 ---#
			# ---- IN THE BELOW LINE 'A2' AND 'B2' ARE RESPECTIVELY THE MAJOR AND MINOR ALLELE OF SNP2 ---#

			# GET HAPLOTYPES FREQUENCIES
			indx <- which(hap.matrix[,1] == 1)
			temp1 <- hap.matrix
			temp1[indx,1] <- 2
			temp2 <-  temp1[,1] + temp1[,2]
			P.AB <- length(which(temp2==0))/num.haps
			P.Ab <- length(which(temp2==2))/num.haps
			P.ab <- length(which(temp2==3))/num.haps
			P.aB <- length(which(temp2==1))/num.haps

			# VERIFY IF HAPLOTYPE FREQUENCIES SUM UP TO 1
			P.AB + P.aB + P.ab + P.Ab

			# GET THE EMPIRICAL ALLELE FREQUENCIES
			P.a <- length(which(hap.matrix[,1] == 1))/num.haps
			P.b <- length(which(hap.matrix[,2] == 1))/num.haps
			P.A <- 1-P.a
			P.B <- 1-P.b

			# MATRIX TO SAVE FREQUENCIES OF THE GENERATED ALLELES
			freqs.matrix <- matrix(0, nrow = 2, ncol = num.snps)
			colnames(freqs.matrix) <- columns
			row.names(freqs.matrix) <- c("MAF", "1-MAF")
			freqs.matrix[1,] <- c(P.a, P.b)
			freqs.matrix[2,] <- c(P.A, P.B)

			# GET THE LEWONTIN D VALUE BETWEEN SNP1 AND SNP2
			D <- P.AB - (P.A * P.B)

			# COMPUTE EMPIRICAL NORMALIZED LEWONTIN D (D')
			if(D < 0){
			   Dprime <- D / min((P.A*P.b), (P.B*P.a))
			}else{
			   Dprime <- D / max(-(P.A*P.B), -(P.a*P.b))
			}

			# COMPUTE THE EMPIRICAL CORRELATION COEFFICIENT R
			r <- D / sqrt(P.A*P.a*P.B*P.b)

			# DISPLAY A SCREEN ALERT AND STOP THE PROCESS IF THE TARGET R IS NOT 
			# COMPATIBLE WITH THE DIFFERENCE IN MAFS


  		snp1.allele <- hap.matrix[,1]
			snp2.allele <- hap.matrix[,2]

			# RETURNED OBJECT
			return(list('alleles'=data.frame(snp1.allele, snp2.allele), 'estimated.r'=r, 'estimated.D'=D, 'estimated.Dprime'=Dprime))
}

