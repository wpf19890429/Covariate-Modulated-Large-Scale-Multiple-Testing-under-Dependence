run_cmlocfdr=function(Pvalue=1,P,X,bases_X=FALSE,K=2,knots=NULL,
	nIter=160,thin=1,burnIn=10,SSA=1,SSG=1,MA=3,MG=3,theoNULL=FALSE,mu,inits=NULL){

	library(pgnorm)
	library(locfdr)
	library(splines)

	## Inputs:
	## P: N vector of P-values or Z scores based on the input from mainfile.R;
	## X: N x Q matrix of covariates
	## nIter: Number of MCMC iterations
	## thin: thinning rate
	## burnIn: burn-in number
	## SSA: increase step-size for alpha draw
	## SS: increase step-size for gamma draw
	## MA: num of multiple try for alpha
	## MG: num of multiple try for gamma
	## mu: origin of gamma distribution, do not change; not used for generalized normal

	if (Pvalue==1){ # input P value
		N=length(P)
		snpid=1:N
		all.complete=complete.cases(cbind(P,X))
		X=X[all.complete,]
		P=P[all.complete]

		P[P==0]=min(P[P>0]);P[P==1]=max(P[P<1])
		Z=-qnorm(P/2);N=length(Z)
		if(!colnames(X)[1]=="Intcpt"){X=cbind(Intcpt=1,X)}
	}else{
		Z=P; #input Z score;
	}
	
	X_bases=NULL
	if(!is.null(bases_X)){
		for(p in 1:length(bases_X)){
			X_p=X[,bases_X[p]]
			range=c(min(X_p),max(X_p))
			delta=(range[2]-range[1])/200
			grid=seq(range[1],range[2]+delta,by=delta)
			knots=c(min(grid),quantile(X_p,.5),quantile(X_p,.67),max(grid))
			phi.mat=bs(grid,knots=knots[2:(length(knots)-1)],degree=3,intercept=FALSE,
				Boundary.knots=c(knots[1],knots[length(knots)]))
			for(k in 1:(K+3)){
				phi.mat[,k]=phi.mat[,k]/(sum(phi.mat[,k])*delta)
			}
			#plot(grid,phi.mat[,1],type="l",ylab="density",xlab="annotation score")
			#title(main="density basis functions")
			#for(k in 2:K){
			#	lines(grid,phi.mat[,k],type="l",col=k+1,lwd=.5)
			#}
			Phi.mat=phi.mat
			for(k in 1:(K+3)){
					for(j in 1:length(grid)){Phi.mat[j,k]=sum(phi.mat[grid<=grid[j],k]*delta)}
			}
			#plot(grid,Phi.mat[,2],type="l")
			#lines(grid,Phi.mat[,3],type="l",col=k)
			X_bases_p=array(NA,dim=c(length(X_p),2))
			for(g in 1:length(grid)){	
				#print(c(p,g))			
				if(length(X_p[abs(X_p-grid[g])<=delta/2])>0){
					tmp=dim(rbind(X_bases_p[abs(X_p-grid[g])<=delta/2,]))[1]
					X_bases_p[abs(X_p-grid[g])<=delta/2,]=cbind(rep(Phi.mat[g,2],tmp),rep(Phi.mat[g,3],tmp))
				}	
			}			
			X_bases=cbind(X_bases,X_bases_p)
		}	
		name=colnames(X)[bases_X[1]]
		colnames(X_bases)=rep(1:2,length(bases_X))
		colnames(X_bases)[1:2]=c(paste(name,"_1",sep=""),paste(name,"_2",sep=""))
		if(length(bases_X)>1){
				for(p in 2:length(bases_X)){
						name=colnames(X)[bases_X[p]]
						colnames(X_bases)[(2*(p-1)+1):(2*(p-1)+2)]=c(paste(name,"_1",sep=""),paste(name,"_2",sep=""))
				}
		}		
		colnames(X_bases)
		X=cbind(Intcpt=1,X_bases,X[,-c(1,bases_X)])
	}
		
	#save(file="data_inputs.R",Z,X,N)

	
	source("cmlFDR_GammaDist.R")
	MCMCfit=cmlFDR_GammaDist(Z,X,nIter=nIter,burnIn=burnIn,thin=thin,SSA=SSA,SSG=SSG,MA=MA,MG=MG,mu=mu,
					theoNULL=theoNULL,inits=inits)

	ALPHA_array=MCMCfit[[1]]
	BETA_array=MCMCfit[[2]]
	GAMMA_array=MCMCfit[[3]]
	SIGMA_SQ_array=MCMCfit[[4]]
	first=1
	Alpha=ALPHA_array[[first]]
	Beta=BETA_array[[first]]
	Gamma=GAMMA_array[[first]]
	if(theoNULL==FALSE){Sigma_sq=SIGMA_SQ_array[[first]]}
	for(j in (first+1):length(ALPHA_array)){
		Alpha=Alpha+ALPHA_array[[j]]
		Beta=Beta+BETA_array[[j]]
		Gamma=Gamma+GAMMA_array[[j]]
		if(theoNULL==FALSE){Sigma_sq=Sigma_sq+SIGMA_SQ_array[[j]]}
	}
	Alpha=cbind(Alpha/length(first:length(ALPHA_array)))
	Beta=Beta/length(first:length(ALPHA_array))
	Gamma=cbind(Gamma/length(first:length(ALPHA_array)))
	if(theoNULL==FALSE){Sigma_sq=Sigma_sq/length(first:length(ALPHA_array))}
	if(theoNULL==TRUE){Sigma_sq=1}
			
	f0<-2*dnorm(abs(Z),mean=0,sd=Sigma_sq^.5)
	f1<-f0 
	for (i in 1:length(Z)){	
		X_i=rbind(X[i,])
		f1[i]<-2*dgamma(abs(Z[i])-.68,shape=exp(X_i%*%Alpha),rate=Beta)
	}
	p0=1
	p1=exp(X%*%Gamma)
	pi0=p0/(p0+p1)
	pi1=1-pi0
	cmfdr=pi0*f0/(pi0*f0+pi1*f1)

#	X_mean=rbind(apply(X,2,mean))
#	f1=2*dpgnorm(Z,p=Beta,mean=0,sigma=exp(X_mean%*%Alpha)*sqrt(gamma(3/Beta)/gamma(1/Beta)))
#	p1=exp(X_mean%*%Gamma)
#	pi0=p0/(p0+p1)
#	pi1=1-pi0
#	fdr=pi0*f0/(pi0*f0+pi1*f1)

#	z.locfdr=locfdr(c(Z,-Z),nulltype=1,plot=0)
#	efron_fdr=z.locfdr$fdr[1:length(Z)]

	results=list()
	results[[1]]=MCMCfit[[1]]
	results[[2]]=MCMCfit[[2]]
	results[[3]]=MCMCfit[[3]]
	results[[4]]=MCMCfit[[4]]	
	results[[5]]=MCMCfit[[5]]	
#	results[[6]]=array(NA,dim=c(N,3))
#	results[[6]][,1]=efron_fdr
#	results[[6]][,2]=fdr
#	results[[6]][,3]=cmfdr
      results[[6]]=cmfdr
				
	return(results)	
}
