##  cmLIS.DP method
cmLIS.DP.func<-function(z_score,X_covariates,niter=300,Rep=6,Gibbs_Iter=4)
{

## USAGE
 # cmLIS.func(z_score,X_covariates,...)

## ARGUMENTS
 # z_score=(z_score[1], ..., z_score[m]): the z-scores
 # X_covariates=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix

## VALUES
 # cmLIS.DP: covariate-modulated local index of signicance (cmLIS)

     z_score<-as.data.frame(z_score)
     MCMC.res<-MCMC_HMM_DP(z_score,X_covariates,nIter=niter,Rep,Gibbs_Iter, theoNULL=FALSE, inits=NULL)
     Indd<-seq(200, niter, 5)
     CmLIS<-NULL 
     for(ii in 1:length(Indd)){
          indss<-Indd[ii]
	    tA<-MCMC.res$A[[indss]]
	    tprob<-MCMC.res$sample_prob[[indss]]
	    tsigma_0<-MCMC.res$sigma_0[[indss]]
	    talpha<-MCMC.res$alpha[[indss]]
	    tsigma_1<-MCMC.res$sigma_1[[indss]]

          bwfw.cmLIS<-bwfw_HMM(z_score[,1], X_covariates, A=tA, prob=tprob, sigma_0=tsigma_0, alpha=talpha, sigma_1=tsigma_1)
          tcmLIS<-bwfw.cmLIS$cmLIS
	    CmLIS<-rbind(CmLIS,tcmLIS)
     }	
     cmLIS<-colMeans(CmLIS)    
     dat<-list(cmLIS.DP=cmLIS)
     dat
}