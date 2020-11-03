#' Pseudo Marginal Bayesian Tensor Clustering for Precision Matrices
#'
#' This function performs clustering of Tensor Normal samples with zero mean
#' and Different Precision Matrices throught the use of Pseudo Data Marginal
#' distributions. Based on the Pseudo-Marginal the procedure performs Chinese
#' Restaurant Stochastic Search based on a Dirichlet Process Prior.
#'
#' @param iP       A list specifying the Initial Partition.
#' @param kPseCov  List of Vecorized k-Matrization Pseudo Covariances for each
#'   dimension and for each observation.
#' @param theta    Dirichlet Process Prior Parameter.
#' @param burnin   Number of sample to burn in in the MCMC.
#' @param nmcmc    Number of MCMC samples desired as output. That is without
#'   considering the burn-in period.
#' @param progress Boolean indicating the use of the progress bar. By default
#'   it is set true.
#'
#' @return A list containing samples from the posterior partition.
#'
#' @author Rene Gutierrez Marquez
#'
#' @export
#'

PMBTC <- function(iP,
                  theta,
                  kPseCov,
                  burnin   = 0,
                  nmcmc    = 100,
                  progress = TRUE){
  ### Number of Observations
  n <- ncol(kCov[[1]])
  ### Partition Samples
  sP <- list()
  ### Partition Initialization
  P <- iP
  ### Progress Bar
  if(progress){
    pb <- txtProgressBar(min     = 0,
                         max     = 1,
                         initial = 0,
                         style   = 3,
                         width   = 72)
  }
  ### For Every Stochastic Search Iteration
  for(s in 1:(burnin + nmcmc)){
    ### For all observations
    for(i in 1:n){
      ## Remove observation for the Partition
      H  <- length(P)
      nP <- P
      ### Looks for the Observation
      for(h in 1:H){
        if(i %in% P[[h]]){
          pos <- h
          break
        }
      }
      ### Removes the Element from the Component
      nP[[pos]] <- nP[[pos]][nP[[pos]] != i]
      ### Checks if the New Component is now empty
      if(length(nP[[pos]]) == 0){
        nP <- nP[-pos]
      }

      ## Computes the Number of Elements in the Removed Partition
      nH <- length(nP)
      nL <- numeric(nH)
      for(h in 1:nH){
        nL[h] <- length(nP[[h]])
      }

      ## Arranges the Data According to the New Partition
      pseCov <- tenCovPar(P = nP, pseCov = kCov)
      ## Computes the log-Marginals
      logMar <- numeric(0)
      for(h in 1:nH){
        ### Test Component
        tesCom <- list()
        for(k in 1:K){
          tesCom[[k]] <- cbind(pseCov[[h]][[k]], kCov[[k]][,i])
        }
        logMar <- c(logMar, tenPrePseMar(tesCom) - tenPrePseMar(pseCov[[h]]))
      }
      logMar <- logMar + log(nL)
      logMar <- c(logMar, log(theta))
      #logMar <- c(logMar, -100)

      ## Computes the Probability of Assignment
      pro <- exp(logMar) / sum(exp(logMar))

      ## Assigns the Element to a New Partition
      nh <- sample(1:(nH + 1), size = 1, prob = pro)
      if(nh == (nH + 1)){
        nP[[nh]] <- i
      } else {
        nP[[nh]] <- c(nP[[nh]], i)
      }

      ## Updates the Partition
      P <- nP
    }
    ### Saves the Partition
    if(s > burnin){
      sP[[s - burnin]] <- P
    }
    ### Progress Bar Display
    if(progress){
      setTxtProgressBar(pb    = pb,
                        value = s / (burnin + nmcmc))
    }
  }
  ### Closes the Progress Bar
  if(progress){
    close(pb)
  }

  # Returns the Samples
  return(sP)
}
