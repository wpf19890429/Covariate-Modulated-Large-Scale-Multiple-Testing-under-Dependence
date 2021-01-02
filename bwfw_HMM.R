bwfw_HMM<-function(z, X, A, prob, sigma_0, alpha, sigma_1)
{

###################################################################

## USAGE
 # bwfw_HMM(z, X, pii, A, prob, ...)

## ARGUMENTS
 # z=(z[1], ..., z[m]): the z-scores
 # X=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix
 # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
 # sigma_0: standard deviation of null density
 # alpha=(alpha[1,1],...,alpha[1,p]\\ ... \\ alpha[L,1],...,alpha[L,p])
 # sigma_1: standard deviation of non-null density


## DETAILS
 # bwfw_HMM calculates values for backward, forward variables, probabilities of hidden states, 
 # --cmLIS values and etc. 

## VALUES
 # alpha: rescaled forward variables
 # beta: rescaled backward variables
 # cmLIS: cmLIS values




################################################################## 

## Initialize
    z<-as.numeric(z)
    X<-as.matrix(X)
    NUM<-length(z)
    L<-length(prob)
    pii<-c(0.5, 0.5)

## Densities
    f0z<-dnorm(z, mean=0, sd=sigma_0)
    f1z<-rep(0, NUM)
    if(L>1)
    {   
         for (i in 1:L)
         {
             f1z<-f1z+prob[i]*dnorm(z, mean=X%*%matrix(alpha[i, ], ncol=1), sd=sigma_1[i])
         }
    }else
    f1z<-f1z+prob*dnorm(z, mean=as.vector(X%*%matrix(alpha[1, ], ncol=1)), sd=sigma_1)
    
## the backward-forward procedure

# a. the forward variables
# --rescaled 
    alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

# scaling variable c_0
    c0<-rep(0, NUM)
    alpha[1, 1]<-pii[1]*f0z[1]
    alpha[1, 2]<-pii[2]*f1z[1]

# rescaling alpha
    c0[1]<-1/sum(alpha[1, ])
    alpha[1, ]<-c0[1]*alpha[1, ]

    for (k in 1:(NUM-1))
    { 
      alpha[k+1, 1]<-(alpha[k, 1]*A[1, 1]+alpha[k, 2]*A[2, 1])*f0z[k+1]
      alpha[k+1, 2]<-(alpha[k, 1]*A[1, 2]+alpha[k, 2]*A[2, 2])*f1z[k+1]
# rescaling alpha
      c0[k+1]<-1/sum(alpha[k+1, ])
      alpha[k+1, ]<-c0[k+1]*alpha[k+1, ]
    }

# b. the forward variables
# --rescaled
    beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

    beta[NUM, 1]<-c0[NUM]
    beta[NUM, 2]<-c0[NUM]

    for (k in (NUM-1):1)
    { 
      beta[k, 1]<-A[1, 1]*f0z[k+1]*beta[k+1, 1]+A[1, 2]*f1z[k+1]*beta[k+1, 2]
      beta[k, 2]<-A[2, 1]*f0z[k+1]*beta[k+1, 1]+A[2, 2]*f1z[k+1]*beta[k+1, 2]
# rescaling beta
# using the same scaling factors as alpha 
      beta[k, ]<-c0[k]*beta[k, ]
    }

# c. cmLIS values
# --original
# --the same formulae hold for the rescaled alpha and beta

    cmLIS<-rep(0, NUM)
    for (k in 1:NUM)
    { 
      q1<-alpha[k, 1]*beta[k, 1]
      q2<-alpha[k, 2]*beta[k, 2]
      cmLIS[k]<-q1/(q1+q2)
    }

# d. return the results of the bwfw proc.

    bwfw.res<-list(bw=alpha, fw=beta, cmLIS=cmLIS, c0=c0)
    return(bwfw.res)
  
}




