rdata_HMM1<-function(NUM, pii, A, sigma_0, alpha, sigma_1, X)
{

## USAGE
 # rdata_HMM1(n, pi, A, ...)

## ARGUEMENTS
 # NUM: the number of multiple hypotheses
 # pii=(pii[1], pii[2]): the initial state distribution 
 # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
 # sigma_0: standard deviation of the null density
 # alpha=(alpha_0, alpha_1,..., alpha_p)^T: parameter vector with
 # mu(x)=x_1*alpha_1+...+x_p*alpha_p
 # sigma_1: standard deviation of the non-null density
 # X=(x[1,1], ..., x[1,p] \\ ... \\ x[n,1], ..., x[n,p]): the design matrix


## VALUES
 # rdata_normal1 generates random variables 
 # from a two-group model via a hidden markov model
 # z: continuous observed data
 # theta: binary unobserved states

## Initialize
    theta<-rep(0, NUM)
    z<-rep(0, NUM)

## Generating states
 # initial state
theta[1]<-rbinom(1, 1, pii[2])
 # other states
    for (i in 2:NUM)
    {
      if (theta[i-1]==0)
         theta[i]<-rbinom(1, 1, A[1, 2])
      else
         theta[i]<-rbinom(1, 1, A[2, 2])
    }

## Generating observations
    for (i in 1:NUM)
    {
      if (theta[i]==0)
      {
         z[i]<-rnorm(1, 0, sigma_0)
      }
      else
      { 
         z[i]<-rnorm(1, X[i, ]%*%alpha, sigma_1)
      }
    }
    data<-list(s=theta, o=z)
    return (data)
}



