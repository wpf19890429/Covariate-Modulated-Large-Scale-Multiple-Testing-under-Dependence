## library or source the needed packages and functions 
library(MASS)  # for mvnorm function
library(coda)
library(MCMCpack)  # for rdirichlet function
library(locfdr)
library(pgnorm)
library(splines)
library(tmvtnorm)
library(mnormt)
library(magic)
library(arm)
library(pscl)

source("rdata_HMM.R")
source("Sample_theta.R")
source("MCMC_HMM.R")
source("bwfw_HMM.R")
source("mt.hmm.R.txt")
source("em.hmm.R.txt")
source("bwfw.hmm.R.txt")
source("mt.bh.R.txt")
source("run_cmlocfdr.R")
source("cmlFDR_GammaDist.R")
source("fcns.R")
source("cmLIS.func.R")
source("cmLIS.DP.func.R")
source("ProgressBar.R") # record the progress for the methods
source("MCMC_HMM_DP.R")

##################################################################################################

