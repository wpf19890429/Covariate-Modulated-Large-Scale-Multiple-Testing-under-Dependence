DP_alpha_sigma2<-function(Z,X,Rep=200,Gibbs_Iter=10,alpha=0.1,m=5,a=1,b=1,s=1)
{

###################################################################

## USAGE
 # DP_alpha_sigma2(Z,X,Rep=200,...)


## ARGUMENTS 
 # Z=(z[1], ..., z[m']): the z-scores 
 # X=(x[1,1], ..., x[1,p] \\ ... \\ x[m',1], ..., x[m',p]): the design matrix
 # Rep: iteration steps of MCMC algorithm
 # alpha: the alpha parameter of Dirichlet Process DP(alpha,G0)
 # m:  the number of auxiliary components
 # a, b: hyperparameters of G0 with shape=a, scale=b
 # Gibbs_Iter: number of iteration of Gibbs algorithm when sampling (alpha, sigma2)


## VALUES
 # prob=(prob[1],..., prob[L]): mixing ratio of each component
 # sample_alpha: L by p matrix
 # sample_sigma1: 1 by L vector


	
	set.seed(s)
	n=length(Z)
	if(sum(X[,1])!=n){
		X<-cbind(1,X)
	}
	p<-ncol(X)-1

	# step0: Initialize 
	Index1=sample(1:n,n/2,replace = FALSE); Index2=setdiff(1:n,Index1)
	cc=rep(0,n);cc[Index1]=1;cc[Index2]=2

	# (p+1)*L matrix
	Beta<-t(mvrnorm(n = 2, mu=rep(0,p+1), Sigma=diag(p+1), tol = 1e-6, empirical = FALSE))
	Sigma2<-rep(1,2)


	
	# MCMC starts
	rep=1
	while(rep<Rep+1){
		#cat("current MCMC rep=",rep,"/",Rep,"\n",sep="")
	
		# step1: drwa c
		for(i in 1:n){
			#cat("Wrong step1","\n")
			K_minus=length(unique(cc[-i]))
			h=K_minus+m
		
			Beta_m<-t(mvrnorm(n = m, mu=rep(0,p+1), Sigma=diag(p+1), tol = 1e-6, empirical = FALSE))
			Sigma2_m<-rinvgamma(m, shape=a, scale = b) #rate=1/scale

			if(K_minus<length(unique(cc))){ #c[i]!=c[j] for any j!=i
				# relabel c[-i], Beta, Sigma2 according to 1:K_minus
				for(c_test in sort(unique(cc))){
					if(c_test>cc[i]){
						cc[cc==c_test]<-c_test-1
					}
				}
				update_ind<-sort(setdiff(1:ncol(Beta),cc[i]))
				Beta_K_minus_Plus_1<-Beta[,cc[i]]
				Sigma2_K_minus_Plus_1<-Sigma2[cc[i]]

				Beta<-Beta[,update_ind]
				Beta<-matrix(Beta,nrow=(p+1))
				Sigma2<-Sigma2[update_ind]
				cc[i]<-K_minus+1
				Beta<-cbind(Beta,Beta_K_minus_Plus_1)
				Sigma2<-c(Sigma2,Sigma2_K_minus_Plus_1)
			
				Beta_m<-Beta_m[,-m]
				Sigma2_m<-Sigma2_m[-m]
			}		
			Beta<-cbind(Beta,Beta_m)
			Sigma2<-c(Sigma2,Sigma2_m)

			# draw c[i] from {1, 2, ..., h}
			prob<-rep(0,h)
			for(c_sample in 1:K_minus){
				n_minusi_c=length(which(cc[-i]==c_sample))
				prob[c_sample]=n_minusi_c/(n-1+alpha)*dnorm(Z[i], mean = as.numeric(X[i,]%*%Beta[,c_sample]), sd = sqrt(Sigma2[c_sample]), log = FALSE)		
			}
			for(c_sample in (K_minus+1):h){
				prob[c_sample]=(alpha/m)/(n-1+alpha)*dnorm(Z[i], mean = as.numeric(X[i,]%*%Beta[,c_sample]), sd = sqrt(Sigma2[c_sample]), log = FALSE)		
			}
			prob<-prob/sum(prob)
			cc[i]<-sample(1:h,1,prob=prob)

			# change the state Beta, Sigma2 to contain only the beta, sigma2 that are now associated with one or more observations
			if(cc[i]<K_minus+1){
				Beta<-Beta[,1:K_minus]; Beta<-matrix(Beta,nrow=(p+1))
				Sigma2<-Sigma2[1:K_minus]
			}else{
				Beta<-Beta[,c(1:K_minus,cc[i])];Beta<-matrix(Beta,nrow=(p+1))
				Sigma2<-Sigma2[c(1:K_minus,cc[i])]
				cc[i]<-K_minus+1
			}
		}
	
		## step2: sample {Beta,Sigma2} by Gibbs method
		for(cc_sample in 1:ncol(Beta)){
			#cat("Wrong step2","\n")
			inds<-which(cc==cc_sample)

			if(length(inds)==1){
				tXX<-t(t(X[inds,]))%*%X[inds,]
				sum_zX<-t(t(Z[inds]*(X[inds,])))
			}else{
				tXX<-t(X[inds,])%*%X[inds,]
				sum_zX<-t(X[inds,])%*%Z[inds]
			}
		

			for(Gibbs_iter in 1:Gibbs_Iter){
				# given sigma2 to sample beta form beta|sigma2,Z[inds]
				SIGMA_Inv<-diag(p+1)+1/(Sigma2[cc_sample]+1e-8)*tXX
				SIGMA<-ginv(SIGMA_Inv)
				Mu<-1/(Sigma2[cc_sample]+1e-8)*(SIGMA%*%sum_zX)
				Beta[,cc_sample]<-t(mvrnorm(n = 1, mu=Mu, Sigma=SIGMA, tol = 1e-6, empirical = FALSE))
			
				# given beta to sample sigma2 form sigma2|beta,Z[inds]
				shape=a+length(inds)/2
		
				temp_s<-0
				for(tem_i in 1:length(inds)){
					temp_s<-temp_s+(Z[inds[tem_i]]-X[inds[tem_i],]%*%Beta[,cc_sample])^2
				}
				scale=as.numeric(b+1/2*temp_s)
				Sigma2[cc_sample]<-rinvgamma(1, shape=shape, scale = scale) #rate=1/scale
			}
		}
		rep=rep+1

	}

	L<-length(unique(cc))
	I_L<-diag(L)
	c_matrix=I_L[cc,]
	prob=colMeans(c_matrix)

	results<-list(cc=cc,Beta=Beta,Sigma2=Sigma2, prob=prob,sample_alpha=t(Beta),sample_sigma1=sqrt(Sigma2))
	return(results)
}


DP_alpha_sigma2_Initial<-function(Z,X,cc,Beta,Sigma2,Rep=200,Gibbs_Iter=10,alpha=0.1,m=5,a=1,b=1,s=1)
{

###################################################################

## USAGE
# DP_alpha_sigma2_Initial(Z,X,cc,Beta,Sigma2,Rep=200,...), where cc,Beta and Sigma2 are as the initial input of this function


## ARGUMENTS 
# Z=(z[1], ..., z[m']): the z-scores 
# X=(x[1,1], ..., x[1,p] \\ ... \\ x[m',1], ..., x[m',p]): the design matrix
# cc=(cc[1], ..., cc[m']): the variable that indicate which group each observation belongs to 
# Beta: L by p matrix
# Sigma2: 1 by L vector
# Rep: iteration steps of MCMC algorithm
# alpha: the alpha parameter of Dirichlet Process DP(alpha,G0)
# m:  the number of auxiliary components
# a, b: hyperparameters of G0 with shape=a, scale=b
# Gibbs_Iter: number of iteration of Gibbs algorithm when sampling (alpha, sigma2)


## VALUES
# prob=(prob[1],..., prob[L]): mixing ratio of each component
# sample_alpha: L by p matrix
# sample_sigma1: 1 by L vector

	
	set.seed(s)
	n=length(Z)
	if(sum(X[,1])!=n){
		X<-cbind(1,X)
	}
	p<-ncol(X)-1

	# step0: Initialize 
	cc<-cc
	Beta<-Beta;Beta<-matrix(Beta,nrow=(p+1))
	Sigma2<-Sigma2


	
	# MCMC starts
	rep=1
	while(rep<Rep+1){
		#cat("current MCMC rep=",rep,"/",Rep,"\n",sep="")
	
		# step1: drwa c
		for(i in 1:n){
			#cat("Wrong step1","\n")
			K_minus=length(unique(cc[-i]))
			h=K_minus+m
		
			Beta_m<-t(mvrnorm(n = m, mu=rep(0,p+1), Sigma=diag(p+1), tol = 1e-6, empirical = FALSE))
			Sigma2_m<-rinvgamma(m, shape=a, scale = b) #rate=1/scale

			if(K_minus<length(unique(cc))){ #c[i]!=c[j] for any j!=i
				# relabel c[-i], Beta, Sigma2 according to 1:K_minus
				for(c_test in sort(unique(cc))){
					if(c_test>cc[i]){
						cc[cc==c_test]<-c_test-1
					}
				}
				update_ind<-sort(setdiff(1:ncol(Beta),cc[i]))
				Beta_K_minus_Plus_1<-Beta[,cc[i]]
				Sigma2_K_minus_Plus_1<-Sigma2[cc[i]]

				Beta<-Beta[,update_ind]
				Beta<-matrix(Beta,nrow=(p+1))

				Sigma2<-Sigma2[update_ind]
				cc[i]<-K_minus+1
			
				Beta<-cbind(Beta,Beta_K_minus_Plus_1)
				Sigma2<-c(Sigma2,Sigma2_K_minus_Plus_1)
			
				Beta_m<-Beta_m[,-m]
				Sigma2_m<-Sigma2_m[-m]
			}		
			Beta<-cbind(Beta,Beta_m)
			Sigma2<-c(Sigma2,Sigma2_m)

			# draw c[i] from {1, 2, ..., h}
			prob<-rep(0,h)
			for(c_sample in 1:K_minus){
				n_minusi_c=length(which(cc[-i]==c_sample))
				prob[c_sample]=n_minusi_c/(n-1+alpha)*dnorm(Z[i], mean = as.numeric(X[i,]%*%Beta[,c_sample]), sd = sqrt(Sigma2[c_sample]), log = FALSE)		
			}
			for(c_sample in (K_minus+1):h){
				prob[c_sample]=(alpha/m)/(n-1+alpha)*dnorm(Z[i], mean = as.numeric(X[i,]%*%Beta[,c_sample]), sd = sqrt(Sigma2[c_sample]), log = FALSE)		
			}
			prob<-prob/sum(prob)
			cc[i]<-sample(1:h,1,prob=prob)

			# change the state Beta, Sigma2 to contain only the beta, sigma2 that are now associated with one or more observations
			if(cc[i]<K_minus+1){
				Beta<-Beta[,1:K_minus];Beta<-matrix(Beta,nrow=(p+1)); Sigma2<-Sigma2[1:K_minus]
			}else{
				Beta<-Beta[,c(1:K_minus,cc[i])];Beta<-matrix(Beta,nrow=(p+1))
				Sigma2<-Sigma2[c(1:K_minus,cc[i])]
				cc[i]<-K_minus+1
			}
		}
	
		## step2: sample {Beta,Sigma2} by Gibbs method
		for(cc_sample in 1:ncol(Beta)){
			#cat("Wrong step2","\n")
			inds<-which(cc==cc_sample)

			if(length(inds)==1){
				tXX<-t(t(X[inds,]))%*%X[inds,]
				sum_zX<-t(t(Z[inds]*(X[inds,])))
			}else{
				tXX<-t(X[inds,])%*%X[inds,]
				sum_zX<-t(X[inds,])%*%Z[inds]
			}
		

			for(Gibbs_iter in 1:Gibbs_Iter){
				# given sigma2 to sample beta form beta|sigma2,Z[inds]
				SIGMA_Inv<-diag(p+1)+1/(Sigma2[cc_sample]+1e-8)*tXX
				SIGMA<-ginv(SIGMA_Inv)
				Mu<-1/(Sigma2[cc_sample]+1e-8)*(SIGMA%*%sum_zX)
				Beta[,cc_sample]<-t(mvrnorm(n = 1, mu=Mu, Sigma=SIGMA, tol = 1e-6, empirical = FALSE))
			
				# given beta to sample sigma2 form sigma2|beta,Z[inds]
				shape=a+length(inds)/2
		
				temp_s<-0
				for(tem_i in 1:length(inds)){
					temp_s<-temp_s+(Z[inds[tem_i]]-X[inds[tem_i],]%*%Beta[,cc_sample])^2
				}
				scale=as.numeric(b+1/2*temp_s)
				Sigma2[cc_sample]<-rinvgamma(1, shape=shape, scale = scale) #rate=1/scale
			}
		}
		rep=rep+1

	}

	L<-length(unique(cc))
	I_L<-diag(L)
	c_matrix=I_L[cc,]
	prob=colMeans(c_matrix)

	results<-list(cc=cc,Beta=Beta,Sigma2=Sigma2, prob=prob,sample_alpha=t(Beta),sample_sigma1=sqrt(Sigma2))
	return(results)
}











##################################################################################################################


