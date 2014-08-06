#' 
#' @title Checks if a matrix is positive definite
#' @description  Checks if any of the eigenvalues of the matrix is smaller than the set tolerance value
#' @param matrix input matrix
#' @param tolerance a constant
#' @return A boolean
#' @keywords internal
#' @author Gaye A.
#' 
is.posdef <-
function(matrix, tolerance=0.000001)
{

  # IF ONE OR MORE EIGEN VALUES < TOLERANCE VALUE, MATRIX IS NOT POSITIVE DEFINITE
  negative.eigens <- length(which(eigen(matrix)$values < tolerance))
  if(negative.eigens > 0){FALSE}else{TRUE}
  
}

