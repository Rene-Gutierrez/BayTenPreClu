#' Pseudo Marginal Distribution of Normal Tensor
#'
#' Evaluates the a Pseudo Marginal Distribution of Normal Tensor Data with zero
#' mean and precision matrices are unknown. Note that this is not the actual
#' marginal or proportional to it.
#'
#' @param C A list of matrices, one for each dimension of the tensor containing
#'   the matricization of the tensor by each dimension.
#'
#' @return The value of the pseudo log-marginal distribution.
#'
#' @author Rene Gutierrez Marquez

#' @export

tenPrePseMar <- function(C){
  ### Gets the Number of Dimensions
  K <- length(C)
  ### Gets the Mean of Each Dimension
  mC <- list()
  for(k in 1:K){
    mC[[k]] <- rowMeans(C[[k]])
  }
  ### Computes the log of the marginal
  logMar <- 0
  for(k in 1:K){
    p      <- length(mC[[k]])
    logMar <- logMar + sum((C[[k]] - mC[[k]])^2) / p
  }
  logMar <- logMar / K
  return(-logMar)
}
