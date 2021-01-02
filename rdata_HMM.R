rdata_HMM<-function(NUM, pii, A, sigma_0, prob, alpha, sigma_1, X)
{

## USAGE
 # rdata_HMM(n, pi, A, ...)

## ARGUEMENTS
 # NUM: the number of multiple hypotheses
 # pii=(pii[1], pii[2]): the initial state distribution 
 # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
 # sigma_0: standard deviation of null density
 # prob: proportion of mixed components under the alternative hypothesis
 # alpha=(alpha_0, alpha_1,..., alpha_p)^T: parameter vector with
 # mu(x)=x_1*alpha_1+...+x_p*alpha_p
 # sigma_1: standard deviation of non-null density
 # X=(x11, ..., x1p \\ ... \\ xn1, ..., xnp): design matrix


## VALUES
 # rdata_normal generates random variables 
 # from a two-group model via a hidden markov model
 # o=z: continuous observed data
 # s=theta: binary unobserved states

## Initialize
    theta<-rep(0, NUM)
    z<-rep(0, NUM)
    L<-length(sigma_1)


## Generating the states
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

## Generating the observations
    for (i in 1:NUM)
    {
      if (theta[i]==0)
      {
         z[i]<-rnorm(1, 0, sigma_0)
      }
      else
      { 
         c<-sample(1:L, 1, prob=prob)
         z[i]<-rnorm(1, mean=X[i, ]%*%alpha[c, ], sd=sigma_1[c])        
      }
    }

    data<-list(s=theta, o=z)
    return(data)

}

