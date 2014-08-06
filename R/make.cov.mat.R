#' 
#' @title Generates the covariance matrix required to achieved the desired LD
#' @description  Finds the covariance values required to achieve the specified frequency 
#' of major allele haplotype
#' @param cor.mat Correlation matrix
#' @param freqs Major allele frequencies of the two snps
#' @param is.add Models of the two variants
#' @return A 2X2 covariance matrix
#' @keywords internal
#' @author Gaye A.
#' @references Montana, G. 2005, \code{HapSim: a simulation tool for 
#' generating haplotype data with pre-specified allele frequencies 
#' and LD coefficients.}, Bioinformatics, \bold{vol. 21 (23)}, 
#' pp.4309-4311.
#' 
make.cov.mat <-
function(cor.mat, freqs){

# COMPUTES CUMULATIVE DISTRIBUTION FUNCTION FOR A BIVARIATE NORMAL DISTRIBUTION
cdfbn <- function(a,b,c)
{
  a <- a
  b <- b
  c <- c
  aa = c(0.325303,0.4211071,0.1334425,0.006374323)
  bb = c(0.1337764,0.6243247,1.3425378,2.2626645)
  PI <- 3.141592654
  res <- 0

  if(c > 0.9999){
    output <- pnorm(min(a,b))
  }
  if(c < -0.9999){
    output = max(0, pnorm(a) - pnorm(-b))
  }
  if(a*b*c <= 0){
     if(a <= 0 & b <= 0 & c <= 0){
       for(i in 1:4){
         for(j in 1:4){
           cc <- a / sqrt(2*(1 -c*c))
           dd <- b / sqrt(2*(1-c*c))
           ee <- exp(cc*(2*bb[i]-cc) + dd*(2*bb[j]-dd) + 2*c*(bb[i]-cc)*(bb[j]-dd))
           res <- res + (aa[i] * aa[j] * ee)
         }
       }
     res <- res * sqrt(1-c*c)/PI
     output <- res
     }
     if(a <= 0 & b*c >= 0){
       for(i in 1:4){
         for(j in 1:4){
           cc <- a / sqrt(2*(1 -(-c)*(-c)))
           dd <- -b / sqrt(2*(1-(-c)*(-c)))
           ee <- exp(cc*(2*bb[i]-cc) + dd*(2*bb[j]-dd) + 2*(-c)*(bb[i]-cc)*(bb[j]-dd))
           res <- res + (aa[i] * aa[j] * ee)
         }
       }
       res <- res * sqrt(1-(-c)*(-c))/PI
       output <- pnorm(a) - res
     }
     if(b <= 0 & c >= 0){
       for(i in 1:4){
         for(j in 1:4){
           cc <- -a / sqrt(2*(1 -(-c)*(-c)))
           dd <- b / sqrt(2*(1-(-c)*(-c)))
           ee <- exp(cc*(2*bb[i]-cc) + dd*(2*bb[j]-dd) + 2*(-c)*(bb[i]-cc)*(bb[j]-dd))
           res <- res + (aa[i] * aa[j] * ee)
         }
       }
       res <- res * sqrt(1-(-c)*(-c))/PI
       output <- pnorm(b) - res
     }else{
       for(i in 1:4){
         for(j in 1:4){
           cc <- -a / sqrt(2*(1 -c*c))
           dd <- -b / sqrt(2*(1-c*c))
           ee <- exp(cc*(2*bb[i]-cc) + dd*(2*bb[j]-dd) + 2*c*(bb[i]-cc)*(bb[j]-dd))
           res <- res + (aa[i] * aa[j] * ee)
         }
       }
       res <- res * sqrt(1-(-c)*(-c))/PI
       output <- pnorm(a) + pnorm(b) - 1 +res
     }
  }else{
     denominator = sqrt(a * a - 2 * c * a * b + b * b)
     prob1 = (c * a - b) * sign(a) / denominator
     prob2 = (c * b - a) * sign(b) / denominator
     output = cdfbn(a, 0, prob1) + cdfbn(b, 0, prob2) - (1 - sign(a) * sign(b)) / 4
  }
  if(output < 0){
    output <- 0
  }
  return (output)
}

# SOLVES THE EQUATION FOR THE CORRELATION: WHAT COVARIANCE VALUE IS EQUAL TO THIS C.D.F?
find.cor <- function(a,b,c)
{
  aa <- -1
  bb <- 1
  precision <- 0.00001
 
  while ((bb-aa) > precision){
		  output <- (aa+bb)/2
		  if(cdfbn(a,b,output) == c){
     break
		  }
		  if(cdfbn(a,b,output) < c){
			  aa <- output
    }
    if(cdfbn(a,b,output) > c){
			  bb <- output
	   }
  }
	return (output)
}

# computes the required covariance matrix
quantiles <- abs(qnorm(freqs))
cov.mat <- matrix(c(0,0,0,0),2,2)
	 for(i in 1:2){
		  for(j in 1:2){
			   aa <- (cor.mat[i,j]*sqrt(freqs[1]*(1-freqs[1])*freqs[2]*(1-freqs[2]))) + (freqs[1]*freqs[2])
   			cov.mat[i,j] <- find.cor(quantiles[1], quantiles[2], aa)
		  }
  }
return (cov.mat)
}

