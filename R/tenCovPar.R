#' Tensor Data Splitter by Partition
#'
#' Partitions
#'
#' @param P      A list containing a partition of indexes.
#' @param pseCov
#'
#' @return A list of lists containing the ...
#'
#' @author Rene Gutierrez Marquez

#' @export

tenCovPar <- function(P, pseCov){
  ### Number of Components
  H <- length(P)
  ### Number of Dimensions
  K <- length(pseCov)
  ### For every Component Creates a List
  Pcov <- list()
  for(h in 1:H){
    Pcov[[h]] <- list()
    for(k in 1:K){
      Pcov[[h]][[k]] <- as.matrix(pseCov[[k]][, P[[h]]])
    }
  }
  ### Returns the Partition
  return(Pcov)
}
