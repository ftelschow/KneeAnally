################################# Simulations Article ##########################
# 
#  This file contains the simulations of the simultaneous covering rate for
#  the article ""
#
################################################################################

args <- commandArgs(trailingOnly = TRUE)
alpha <- as.numeric(args[[1]])
i <- as.numeric(args[[2]])

library(MASS)
library(expm)
library(KneeMotionAnalytics)

# Change to your path you want to save the simulations
path = "/home/drtea/Research/Rpackages/KneeAnally/Reproduce_Simulations_Figures_article"

M = 2
noise <- c("toyNoise1", "toyNoise2", "toyNoise3")

simulationResults <- data.frame()

# for(noise in c("toyNoise1", "toyNoise2", "toyNoise3")){
  for(sigma in c(0.005, 0.05, 0.1, 0.6)){
    sigma1 <- sigma2 <- sigma3 <- sigma
    for(sName in c("const", "sin")){
      for(corr in c(FALSE, TRUE)){
        for(nSamp in c(10, 15, 30)){
          simulationResults <- rbind(simulationResults,
                                     CovRate.ExpModelSim(
                                       M               = M,
                                       nSamp           = nSamp,
                                       alpha           = alpha,
                                       noise           = noise[i],
                                       sigma1, sigma2, sigma3,
                                       sigmaName       = sName,
                                       timeGrid        = seq(0,1,length.out=100),
                                       corr            = corr
                                     )
                                     )
        }
      }
    }
  }
# }

save(simulationResults, file=paste("path","SO3CovRateResults",100*alpha, noise[i],".RData", sep=""))