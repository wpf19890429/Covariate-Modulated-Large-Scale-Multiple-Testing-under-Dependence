Sample_theta1<-function(z, X, pii, A, sigma_0, alpha, sigma_1)
{

###################################################################

## USAGE
 # Sample_theta1(z, X, pii, A, sigma_0, alpha, sigma_1)

## ARGUMENTS
 # z=(z[1], ..., z[m]): the z-scores
 # X=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix
 # pii=(pii[0], pii[1]): the initial state distribution
 # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
 # sigma_0: standard deviation of the null density
 # alpha=(alpha_0, alpha_1,..., alpha_p)^T: parameter vector with
 # mu(x)=x_1*alpha_1+...+x_p*alpha_p
 # sigma_1: standard deviation of the non-null density

## VALUES
 # theta: the sampling of states from posterior distributions
 # S1: S1[t, 1]=Pr(theta[t]=0|z[1],...,z[t],X)
 # S2: S2[t, 1]=Pr(theta[t]=0|z[1],...,z[t-1],X)

## Initialize
    NUM<-length(z)
    q<-ncol(X);p<-q

## Densities
    f0z<-dnorm(z, mean=0, sd=sigma_0)
    f1z<-dnorm(z, mean=X%*%alpha, sd=sigma_1)

# a. S1 and S2
    S1<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
    S2<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

    S2[1, 1]<-pii[1]
    S2[1, 2]<-pii[2]
    S1[1, 1]<-pii[1]*f0z[1]/(pii[1]*f0z[1]+pii[2]*f1z[1])
    S1[1, 2]<-1-S1[1, 1]

    for(i in 2:NUM)
    {
        S2[i, 1]<-A[1, 1]*S1[i-1, 1]+A[2, 1]*S1[i-1, 2]
        S2[i, 2]<-A[1, 2]*S1[i-1, 1]+A[2, 2]*S1[i-1, 2]
        S1[i, 1]<-f0z[i]*S2[i, 1]/(f0z[i]*S2[i, 1]+f1z[i]*S2[i, 2])
        S1[i, 2]<-1-S1[i, 1]
    }

# b. Sampling theta from theta[NUM] to theta[1]
# In the i-th iteration, Poterior_prob=Pr(theta[i]=1|theta[i+1],z[1],...,z[i],X)  
    theta<-rep(0, NUM)
    theta[NUM]<-rbinom(1, 1, S1[NUM, 2])
    for(i in (NUM-1):1)
    {
         q1<-A[2, theta[i+1]+1]*S1[i, 2] 
         q2<-A[1, theta[i+1]+1]*S1[i, 1] 
         Poterior_prob<-q1/(q1+q2)
         theta[i]<-rbinom(1, 1, Poterior_prob)
    }

    data<-list(theta=theta)
    return(data)
}


