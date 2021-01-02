MCMC_HMM=function(z_score, X_covariates, L, nIter=200, theoNULL=FALSE, inits=NULL)
{

  set.seed(1)
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
    z_score<-as.data.frame(z_score)
    z <- z_score[,1]; 
    X <- X_covariates
    z<-as.numeric(z)
    X<-as.matrix(X)
    NUM<-length(z)
    q<-ncol(X);p<-q
    sample_alpha<-array(rep(0, nIter*L*q), dim=c(nIter, L, q))
    sample_sigma_0<-rep(0, nIter)
    sample_sigma_1<-array(0, dim=c(nIter, L))
    sample_prob<-array(0, dim=c(nIter, L))
    sample_A<-list()




## Prior hyperparameters 
    a0_A<-3
    b0_A<-3
    a0<-2
    b0<-1
    sigma_mu<-10



## Parameters initialization
    sample_alpha[1, , ]<-matrix(rep(1, L*q), L, q)
    sample_sigma_0[1]<-1
    sample_sigma_1[1, ]<-rep(1, L)
    pii<-c(0, 1)
    A<-matrix(rep(0.5, 4), 2, 2, byrow=TRUE)
    prob<-rep(1/L, L)
    theta<-Sample_theta(z, X, pii, A, sample_sigma_0[1], prob, as.matrix(sample_alpha[1, , ]), sample_sigma_1[1, ])$theta
    theta1<-rep(0, NUM)
    NUM_theta1<-length(which(theta==1))
    for(j in which(theta==1))
    {     
         theta1[j]<-sample(1:L, 1, prob=rep(1/L, L))
    } 



## Iterations
    pb<-PB(Methods=paste("MCMC_HMM_L=",L,sep=""),Rep=nIter)
    for(i in 2:nIter)
    {
	    pb$tick()
          #print(paste("iteration", i))

 ## Sample sigma_0
          Z0<-z[theta==0]
          sample_sigma_0[i]<-sqrt(1/rgamma(1, shape=sum(1-theta)/2+a0, rate=t(Z0)%*%Z0/2+b0))
          #print(paste("sample_sigma_0[i]", sample_sigma_0[i]))


 ## Sample prob
          m<-rep(0, L)
          for(j in 1:L)
          {
               m[j]<-length(which(theta1==j))
          }
          prob<-rdirichlet(1, rep(3, L)+m)
          sample_prob[i, ]<-prob  

 ## Sample sigma_1 and Sample alpha

          if(L>1)
          {
               for(j in 1:L)
               {
                   Xj<-X[theta1==j, ]
                   Zj<-z[theta1==j]
                   sample_sigma_1[i, j]<-sqrt(1/rgamma(1, shape=(m[j]-q)/2+a0, rate=t(Zj-Xj%*%sample_alpha[i-1, j, ])%*%(Zj-Xj%*%sample_alpha[i-1, j, ])/2+b0))
                   XjtXj_inv<-solve(t(Xj)%*%Xj)
                   alpha_hat<-XjtXj_inv%*%t(Xj)%*%Zj    
                   sample_alpha[i, j, ]<-mvrnorm(1, alpha_hat, sample_sigma_1[i, j]^2*XjtXj_inv)

               }
          }else
          {   
               X1<-X[theta1==1, ]
               Z1<-z[theta1==1]
               sample_sigma_1[i, 1]<-sqrt(1/rgamma(1, shape=(m[1]-q)/2+a0, rate=t(Z1-X1%*%matrix(sample_alpha[i-1, 1, ], ncol=1))%*%(Z1-X1%*%matrix(sample_alpha[i-1, 1, ], ncol=1))/2+b0))
               X1tX1_inv<-solve(t(X1)%*%X1)
               alpha_hat<-X1tX1_inv%*%t(X1)%*%Z1    
               sample_alpha[i, 1, ]<-mvrnorm(1, alpha_hat, sample_sigma_1[i, 1]^2*X1tX1_inv)     
               #print("sample_alpha[i, ]")
               #print(sample_alpha[i, ])    
               #print(paste("sample_sigma_1[i]", sample_sigma_1[i, ]))
          } 
                


 ## Sample theta and theta1
          theta<-Sample_theta(z, X, pii, A, sample_sigma_0[i], prob, as.matrix(sample_alpha[i, , ]), sample_sigma_1[i, ])$theta
          prob1<-array(0, dim=c(L, length(which(theta==1))))
          X1<-X[theta==1, ]
          Z1<-z[theta==1]
          if(L>1)
          {
                for(j in 1:L)
                {
                      prob1[j, ]<-prob[j]*dnorm(Z1, mean=X1%*%matrix(sample_alpha[i, j, ], ncol=1), sd=sqrt(sample_sigma_1[i, j]))   
                }
          }else
          {
                prob1[1, ]<-prob[1]*dnorm(Z1, mean=X1%*%as.vector(sample_alpha[i, 1, ]), sd=sample_sigma_1[i, 1])
          }

          sum_prob<-apply(prob1, 2, sum)
          theta1<-rep(0, NUM)
          k<-1 
          for(j in which(theta==1))
          {          
                theta1[j]<-sample(1:L, 1, replace=TRUE, prob=prob1[, k]/sum_prob[k])
                k=k+1
          } 


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

