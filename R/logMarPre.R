#' Log-Marginal Distribution of a Normal with Unknown Presicion
#'
#' This function computes the log-marginal distribution of multivariate normal
#' data with zero mean and unknown variance.
#'
#' @param y A list containing the observations of the multivariate normal.
#' @param D Prior scale matrix for the Wishart Distribution, note that this is
#'   not inverted inside the trace.
#' @param b Total degrees of freedom. Again, this differs to the traditional
#'   Wishart representation.
#'
#' @return The value of the log-marginal distribution.
#'
#' @author Rene Gutierrez Marquez

#' @export

logMarPre <- function(y, D, b){
  ### Computes the Length
  n <- length(y)
  k <- length(y[[1]])
  ### Computes S
  S <- 0
  for(i in 1:n){
    S <- S + y[[i]] %*% t(y[[i]])
  }
  S <- S / n
  ### Computes the log Marginal
  logMar <- -(n + b + k + 1) * log(det(S + D/n)) / 2
  ### Returns the marginal
  return(logMar)
}
