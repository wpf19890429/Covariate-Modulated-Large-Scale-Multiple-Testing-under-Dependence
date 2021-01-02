MCMC_HMM1=function(z, X, nIter=200, theoNULL=FALSE, inits=NULL)
{

###################################################################

## USAGE
 # MCMC_HMM1(z, X, nIter=200,...)

## ARGUMENTS
 # z=(z[1], ..., z[m]): the z-scores
 # X=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix

## VALUES
 # sigma_0: 1 by nIter vector, sampled standard deviations under null hypothesis 
 # sigma_1: 1 by nIter vector, sampled standard deviations under alternative hypothesis 
 # alpha: nIter*q array 
 # A: nIter*2*2, sampled transition matrix


## Initialize
    NUM<-length(z)
    q<-ncol(X)
    sample_alpha<-array(rep(0, nIter*q), dim=c(nIter, q))
    sample_sigma_0<-rep(0, nIter)
    sample_sigma_1<-rep(0, nIter)
    sample_A<-list()

## Prior hyperparameters 
    a0_A<-3
    b0_A<-3
    a0<-2
    b0<-1

## Parameters initialization
    sample_alpha[1, ]<-rep(0.5, q)
    sample_sigma_0[1]<-1
    sample_sigma_1[1]<-1
    pii<-c(0, 1)
    A<-matrix(rep(0.5, 4), 2, 2, byrow=TRUE)
    theta<-Sample_theta1(z, X, pii, A, sample_sigma_0[1], sample_alpha[1, ], sample_sigma_1[1])$theta


## Iterations
    for(i in 2:nIter)
    {
          #print(paste("iteration", i))

 ## Sample sigma_0
          Z0<-z[theta==0]
          sample_sigma_0[i]<-sqrt(1/rgamma(1, shape=sum(1-theta)/2+a0, rate=t(Z0)%*%Z0/2+b0))
	    #print(paste("sample_sigma_0[i]", sample_sigma_0[i]))

 ## Sample sigma_1
          X1<-X[theta==1, ]
          Z1<-z[theta==1]
          sample_sigma_1[i]<-sqrt(1/rgamma(1, shape=(sum(theta)-q)/2+a0, rate=t(Z1-X1%*%sample_alpha[i-1, ])%*%(Z1-X1%*%sample_alpha[i-1, ])/2+b0))
	    #print(paste("sample_sigma_1[i]", sample_sigma_1[i]))

 ## Sample alpha
          X1tX1_inv<-solve(t(X1)%*%X1)
          alpha_hat<-X1tX1_inv%*%t(X1)%*%Z1
          sample_alpha[i, ]<-mvrnorm(1, alpha_hat, sample_sigma_1[i]^2*X1tX1_inv)
	    #print("sample_alpha[i, ]")
          #print(sample_alpha[i, ])

 ## Sample theta
          theta<-Sample_theta1(z, X, pii, A, sample_sigma_0[i], sample_alpha[i, ], sample_sigma_1[i])$theta

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

    data<-list(sigma_0=sample_sigma_0, sigma_1=sample_sigma_1, alpha=sample_alpha, A=sample_A)
    return(data)
}

