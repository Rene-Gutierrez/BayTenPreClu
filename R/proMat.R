#' Observations Probability Matrix
#'
#' This function computes the probability matrix following Lau and Green (2007)
#'   according to a posterior sample Partition. The posterior probability
#'   matrix that specifies the posterior probablility that 2 observations
#'   belong to the same cluster.
#'
#' @param sP A list containing posterior sample observations.
#' @param n  Number of observations.
#'
#' @return The Probability Matrix.
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

proMat <- function(sP, n){
  # Set-Up
  ## Gets the Number of Samples
  numSam <- length(sP)
  ## Initializes the Probability Marix
  pM <- matrix(data = 0, nrow = n, ncol = n)

  # Computes the Matrix
  ## For each Sample
  for(i in 1:numSam){
    temMat <- matrix(data = 0, nrow = n, ncol = n)
    numPar <- length(sP[[i]])
    for(j in 1:numPar){
      temMat[sP[[i]][[j]], sP[[i]][[j]]] <- 1
    }
    pM <- pM + temMat
  }

  # Returns the Probability Matrix
  return(pM)
}
