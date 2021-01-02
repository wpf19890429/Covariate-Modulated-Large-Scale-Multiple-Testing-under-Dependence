MCMC_HMM_DP=function(z, X, nIter=200,Rep=20,Gibbs_Iter=10, theoNULL=FALSE, inits=NULL)
{
 source("DP.R")
 source("ProgressBar.R")
 #set.seed(1)
###################################################################

## USAGE
 # MCMC_HMM(z, X, nIter=200,...)

## ARGUMENTS
 # z=(z[1], ..., z[m]): the z-scores
 # X=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix
 # L: the number of mixture components

## VALUES
 # sigma_0: 1 by nIter vector, sampled standard deviations under null hypothesis 
 # sample_prob: nIter by L matrix, sampled mixing ratio of the L mixed components
 # sigma_1: nIter by L matrix, sampled standard deviations of the L mixed components under alternative hypothesis 
 # alpha: nIter*L*q array 
 # A: nIter*2*2, sampled transition matrix

## Initialize
    z<-as.data.frame(z)
    z<-as.numeric(z[,1])
    X<-as.matrix(X)
    NUM<-length(z)
    n<-nrow(X);q<-ncol(X)
    sample_alpha<-list()
    sample_sigma_0<-rep(0, nIter)
    sample_sigma_1<-list()
    sample_prob<-list()
    sample_A<-list()
	
## Prior hyperparameters 
    a0_A<-3
    b0_A<-3
    a0<-2
    b0<-1
    sigma_mu<-10

## Parameters initialization
    L<-1
    q<-ncol(X)
    sample_alpha_initial<-matrix(rep(1, L*q), L, q)
    sample_sigma_0[1]<-1
    sample_sigma_1_initial<-rep(1, L)
    pii<-c(0, 1)
    A<-matrix(rep(0.5, 4), 2, 2, byrow=TRUE)
    prob<-rep(1/L, L)
    theta<-Sample_theta(z, X, pii, A, sample_sigma_0[1], prob, sample_alpha_initial, sample_sigma_1_initial)$theta




## Iterations
    perf<-floor(nIter/3)
    #cat("current MCMC_HMM_DP is running...","\n")

    pb<-PB(Methods="MCMC_HMM_DP",Rep=nIter)
    for(i in 1:nIter)
    {
	      pb$tick()
		#if(i%%(perf)==0){cat("current iter=",i,"/",nIter,"\n",sep="")}# show the progress
          

 ## Sample sigma_0
          Z0<-z[theta==0]
          sample_sigma_0[i]<-sqrt(1/rgamma(1, shape=sum(1-theta)/2+a0, rate=t(Z0)%*%Z0/2+b0))
          #print(paste("sample_sigma_0[i]", sample_sigma_0[i]))


 ## Sample prob 
 ## Sample sigma_1 and Sample alpha



          X1<-matrix(X[theta==1, ],ncol=q)
          Z1<-z[theta==1]
	    
	    if(i==1){
              DP.res<-DP_alpha_sigma2(Z1, X1,Rep,Gibbs_Iter,alpha=0.1,m=3,a=1,b=1,s=2)
		  Inds_theta1=which(theta==1)
		  Cmix<-rep(0,n)
		  Cmix[Inds_theta1]<-DP.res$cc

		  Beta0<-t(mvrnorm(n = 1, mu=rep(0,q), Sigma=diag(q), tol = 1e-6, empirical = FALSE))
		  Sigma2_0<-1				
	    }else{
		  Inds_theta1=which(theta==1)
		  cmix=Cmix[Inds_theta1]   
		  Wind0<-which(cmix==0)
     	        if(length(Wind0)>0){
      		WT<-table(cmix[-Wind0])
      		class_name<-names(WT); 
      		class_prob<-as.numeric(WT)/sum(as.numeric(WT))
      		sample_class<-sample(class_name, length(Wind0), prob = class_prob, replace = T)
      		cmix[Wind0]<-as.numeric(sample_class)
		  }
	  
		  cuq0<-sort(setdiff(unique(cmix),0))
		  cuq<-sort(setdiff(unique(cmix),0))

		  BetaT<-DP.res$Beta[,cuq];BetaT<-matrix(BetaT,nrow=q)
		  Sigma2T<-DP.res$Sigma2[cuq]
		  cuq<-1:length(cuq)# relabel cuq
		  cinds<-sort(which(sort(setdiff(unique(cmix),0))-cuq!=0))
		  for(ss in cinds){
			cmix[which(cmix==cuq0[ss])]<-ss
		  }
	
		  if(0%in%cmix){
			LL<-length(Sigma2T)+1
			cmix[cmix==0]<-LL
			Beta<-cbind(BetaT,t(Beta0));Beta<-matrix(Beta,nrow=q)
			Sigma2<-c(Sigma2T,Sigma2_0)
		  }else{
			Beta<-BetaT
			Sigma2<-Sigma2T
			#print("else")
		  }

		  #if(ncol(Beta)!=length(Sigma2)){cat("wrong here","\n",sep="");1:2%*%matrix(0,5,6)}
		  DP.res<-DP_alpha_sigma2_Initial(Z1,X1,cc=cmix,Beta,Sigma2,Rep,Gibbs_Iter,alpha=0.1,m=3,a=1,b=1,s=2)	  
		  Cmix<-rep(0,n)
		  Cmix[Inds_theta1]<-DP.res$cc
	    }
          sample_prob[[i]]<-DP.res$prob
          sample_alpha[[i]]<-DP.res$sample_alpha
          sample_sigma_1[[i]]<-DP.res$sample_sigma1


 

 ## Sample theta and theta1
          theta<-Sample_theta(z, X, pii, A, sample_sigma_0[i], sample_prob[[i]], sample_alpha[[i]], sample_sigma_1[[i]])$theta



 ## Sample transition matrix A
          C_00<-0
          C_01<-0
          C_10<-0
          C_11<-0
          for(k in 1:(NUM-1))
          {
               if(theta[k]+theta[k+1]==2)
                      C_11<-C_11+1
               else if(theta[k]+theta[k+1]==0)
                      C_00<-C_00+1
               else if(theta[k]>theta[k+1])
                      C_10<-C_10+1
          }
          C_01<-NUM-(C_00+C_10+C_11) 
          a00_new<-rbeta(1, shape1=a0_A+C_00, shape2=b0_A+C_01)
          a11_new<-rbeta(1, shape1=a0_A+C_11, shape2=b0_A+C_10)  
          A<-matrix(c(a00_new, 1-a00_new, 1-a11_new, a11_new), 2, 2, byrow=TRUE) 
          sample_A[[i]]<-A
          #print("A")
          #print(A)

 
    }

    data<-list(sigma_0=sample_sigma_0, sample_prob=sample_prob, sigma_1=sample_sigma_1, alpha=sample_alpha, A=sample_A)
    return(data)
}

