#############################################################################

##                        ---- simulation demo ----

#############################################################################

rm(list=ls())
source("rdata_HMM.R")
source("cmLIS.func.R")
source("cmLIS.DP.func.R")
source("mt.hmm.R.txt")
source("library_and_source.R")


## the number of observed z-values
NUM<-5000
## the initial state distribution
pii<-c(0, 1)
## initialize the transition matrx
A<-matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2, byrow=TRUE)

## initialize parameter set of the null and non-null distributions
sigma_0<-1
sigma_1<-c(1, 1)
alpha<-matrix(c(c(0.5, -2, 1), c(1, 2, 3)), 2, 3, byrow=TRUE)
L<-nrow(alpha)
p<-ncol(alpha)
X<-matrix(0, nrow=NUM, ncol=p)
X[, 1]<-rep(1, NUM)
X[, 2]<-rnorm(NUM, 0, 1)
X[, 3]<-rnorm(NUM, 0, 1)
prob<-c(0.8, 0.2)


## Generating the observed z-valuse and the states of hypotheses that are based on
## covariate-modulated hidden Markov models. 
rdata<-rdata_HMM(NUM, pii, A, sigma_0, prob, alpha, sigma_1, X)
z<-rdata$o
theta<-rdata$s

## Calculating the cmLIS multiple testing statistics by using Bayesian MCMC algorithm
res<-cmLIS.func(z, X, L=2,niter=300)

## Calculating the cmLIS.DP multiple testing statistics by using the non-parametric 
## Bayesian MCMC algorithm
res.DP <-cmLIS.DP.func(z, X, niter=300, Rep=6, Gibbs_Iter=4)



## Conducting cmLIS and cmLIS.DP procedures given the pre-specified level is 0.1. 
level <- 0.1
res_cmLIS <- mt.hmm(res$cmLIS, level)$de # reject or accept
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

#######################################################################################
