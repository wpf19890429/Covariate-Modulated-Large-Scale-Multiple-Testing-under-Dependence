##  cmLIS method
cmLIS.func<-function(z_score,X_covariates,L=2,niter=300)
{
## USAGE
 # cmLIS.func(z_score,X_covariates, L=2,...)

## ARGUMENTS
 # z_score=(z_score[1], ..., z_score[m]): the z-scores
 # X_covariates=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix
 # L: the number of mixture components

## VALUES
 # cmLIS: covariate-modulated local index of signicance (cmLIS)
 # BIC.value: BIC value
     z_score<-as.data.frame(z_score)
     p<-ncol(X_covariates)
     m<-nrow(X_covariates)
     try_num=0
     failure <- try("a" + "b",  silent = T)
     while('try-error' %in% class(failure)){
         try_num <- try_num+1
	   failure <- try(MCMC.res<-MCMC_HMM(z_score, X_covariates, L, nIter=niter, theoNULL=FALSE, inits=NULL),  silent = F)
     }
     

     MCMC.A<-matrix(0, 2, 2)
     MCMC.prob<-rep(0, L)
     MCMC.sigma_0<-0
     MCMC.alpha<-matrix(0, nrow=L, ncol=p)
     MCMC.sigma_1<-rep(0, L)
     for(index in seq(200, niter, 5))
     {
           MCMC.A<-MCMC.A+MCMC.res$A[[index]]/length(seq(200, niter, 5))
           MCMC.prob<-MCMC.prob+MCMC.res$sample_prob[index, ]/length(seq(200, niter, 5))
           MCMC.sigma_0<-MCMC.sigma_0+MCMC.res$sigma_0[index]/length(seq(200, niter, 5))
           MCMC.alpha<-MCMC.alpha+MCMC.res$alpha[index, , ]/length(seq(200, niter, 5)) 
           MCMC.sigma_1<-MCMC.sigma_1+MCMC.res$sigma_1[index, ]/length(seq(200, niter, 5))
     }
     bwfw.cmLIS<-bwfw_HMM(z_score[,1], X_covariates, A=MCMC.A, prob=MCMC.prob, sigma_0=MCMC.sigma_0, alpha=MCMC.alpha, sigma_1=MCMC.sigma_1)
     BIC.cmLIS.value<-2*sum(log(bwfw.cmLIS$c0))+(L*((p-1)+3)+1)*log(m)

     ## compute the BIC-value
     cmLIS<-bwfw.cmLIS$cmLIS
     dat<-list(cmLIS=cmLIS,BIC.value=BIC.cmLIS.value) 
     dat    
}