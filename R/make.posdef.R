#' 
#' @title Turns a matrix into a positive definite one
#' @description  Computes the nearest positive definite matrix of a real symmetric matrix
#' @param matrix input matrix
#' @param tolerance a constant
#' @return A positive-definite matrix
#' @keywords internal
#' @author Gaye A.
#' @references N.J. Higham, 1988 \code{Computing a nearest symmetric positive 
#' semidefinite matrix}, Linear Algebra Appl. \bold{vol. 103}, pp.103 118
#' 
make.posdef <-
function(matrix, tolerance=0.000001)
{
    # IF THE INPUT MATRIX IS NOT SQUARED STOP THE PROCESS 
    if(dim(matrix)[1] != dim(matrix)[2]){
      cat("\n Your matrix is not square! :\n")
      stop()
    }

    # COMPUTE THE EIGEN VALUES
    eig <- eigen(matrix)
    eigvals <- eig$values
 
    # COMPUTE THE NEAREST POSITIVE DEFINITE OF MATRIX, USING THE ALGORITHM OF NJ HIGHAM (1988)
    D = 2 * tolerance
    para.max = pmax(0, D - eigvals)
    pos.def.mat = eig$vectors %*% diag(para.max, dim(matrix)[1]) %*% t(eig$vectors)
    return(matrix + pos.def.mat)
}

