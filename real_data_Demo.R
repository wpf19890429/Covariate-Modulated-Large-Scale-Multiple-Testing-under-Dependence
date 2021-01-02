##################################################################################################


## ------------------------------- cmLIS and cmLISDP demo ------------------------------------- ##


##################################################################################################
rm(list = ls())
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


## --------------------------------- load data ----------------------------- ##

## true L=2
z_score<-read.csv("z_score.csv",head=T)
X_covariates<-read.csv("X_covariates.csv",head=T) 
theta<-read.csv("true_state.csv",head=T)$theta
m<-length(theta)


## --------------------------------- run cmLIS ----------------------------- ##

res1<-cmLIS.func(z_score,X_covariates,L=1)
res2<-cmLIS.func(z_score,X_covariates,L=2)
res3<-cmLIS.func(z_score,X_covariates,L=3)

BIC_values<-c(res1$BIC.value,res2$BIC.value,res3$BIC.value)
print(BIC_values)


## --------------------------------- run cmLIS.DP --------------------------- ##

res.DP<-cmLIS.DP.func(z_score,X_covariates,niter=300,Rep=6,Gibbs_Iter=4)

## ------------------------------ summary the results ------------------------ ##

level <- 0.1
res_cmLIS <- mt.hmm(res2$cmLIS, level)$de # reject or accept
res_cmLIS.DP <- mt.hmm(res.DP$cmLIS.DP, level)$de # reject or accept


## compute FDR, FNR and ATP for cmLIS

N10<-length(which(res_cmLIS-theta==1))
N01<-length(which(theta-res_cmLIS==1))
R<-length(which(res_cmLIS==1))+0.0001
S<-m-R
FDR_cmLIS<-N10/R
FNR_cmLIS<-N01/S
ATP_cmLIS<-R-N10


## compute FDR, FNR and ATP for cmLIS.DP

N10<-length(which(res_cmLIS.DP-theta==1))
N01<-length(which(theta-res_cmLIS.DP==1))
R<-length(which(res_cmLIS.DP==1))+0.0001
S<-m-R
FDR_cmLIS.DP<-N10/R
FNR_cmLIS.DP<-N01/S
ATP_cmLIS.DP<-R-N10


## print the results
cat("\n",rep("-",30),"summary the results",rep("-",30),"\n",sep="")

cat("FDR, FNR and ATP of cmLIS are ", FDR_cmLIS," ", FNR_cmLIS," ", ATP_cmLIS,".","\n",sep="")
cat("FDR, FNR and ATP of cmLIS.DP are ", FDR_cmLIS.DP," ", FNR_cmLIS.DP," ", ATP_cmLIS.DP,".","\n",sep="")

cat("\n",rep("-",30),"end",rep("-",30),"\n",sep="")
##################################################################################################
