#######################
#functions (fcns.R);
#######################

log_p_alpha=function(Alpha,Z,W,B,SA,MU){
	gammafcn=lgamma(exp(W%*%Alpha));
	#if (gammafcn == Inf || is.na (gammafcn) ) {print(paste("Gamma function",gammafcn));}
	log_p=sum(exp(W%*%Alpha)*log(Z-MU)-gammafcn)+sum(exp(W%*%Alpha))*log(B)-0.5*t(Alpha)%*%solve(SA)%*%Alpha
	#print(log_p)
	return(log_p)
}


## Draw Alpha,log-scale MH alg:multiple-try;

Draw_Alpha_log_M_mu=function(Alpha,Z,X,Beta,Phi,df,Multiple,Sigma_Alpha,SSA,mu)
{					
			
		log_p_Alpha_star=rep(0,Multiple)
		log_p_Alpha_2star=rep(0,Multiple)
		p=rep(0,Multiple)
		den=0;num=0;

		sigma=solve(t(X)%*%X);
		diag(sigma)=diag(sigma)*SSA;

		#if(det(sigma) == 0) {print(sigma); break}
		
		logp=log_p_alpha(Alpha,Z,X,B=Beta,SA=Sigma_Alpha,MU=mu);
		if( logp == Inf | logp == -Inf | is.na(logp)) {alpha=Alpha}
		else {alpha <-optim(Alpha,Z=Z,W=X,B=Beta,SA=Sigma_Alpha,MU=mu,log_p_alpha,method="Nelder-Mead",
				hessian=FALSE,control=list(maxit=10,fnscale=-1))$par}

		Alpha_star=t(rmvt(n=Multiple,alpha,sigma=sigma,df=df))
			
			for (i in 1:Multiple){
				log_p_Alpha_star[i]=log_p_alpha(Alpha_star[,i],Z,X,Beta,Sigma_Alpha,mu)
			}
			
			#control overfloat, -max(log_p_Alpha_star);
			p=exp(log_p_Alpha_star-max(log_p_Alpha_star))/sum(exp(log_p_Alpha_star - max(log_p_Alpha_star)))
			
			#in case there is still overfloat;
			p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
			j=sample(c(1:Multiple),1,prob=p);							
			
			Alpha_2star=t(rmvt(n=Multiple-1,Alpha_star[,j],sigma=sigma,df=df))
			Alpha_2star <-cbind(Alpha_2star,Alpha)
			
			for (i in 1:Multiple){
				log_p_Alpha_2star[i]=log_p_alpha(Alpha_2star[,i],Z,X,Beta,Sigma_Alpha,mu)
			}
			
			#control overfloat
			num=sum(exp(log_p_Alpha_star -max(log_p_Alpha_star)))
			den=sum(exp(log_p_Alpha_2star -max(log_p_Alpha_star)))
											
			rho=min(1,num/den)
			
			#in case overfloat again
			if(is.na(rho)) {rho=0.5};
			
			accp=0;
			u=runif(1)
			if(u<rho){
				Alpha=Alpha_star[,j]
				accp=1;
			}
		
		
		return(list(par=Alpha, accp=accp))
}

log_p_gamma=function(Gamma,Z,X,SG,phi){
	log_p=sum(phi*(X%*%Gamma)-log(1+exp(X%*%Gamma)))-0.5*t(Gamma)%*%solve(SG)%*%Gamma
	return(log_p)
}

##Draw Gamma, log scale, multiple-try MH;
Draw_Gamma_log_M=function(Gamma,Z,X,Phi,df,Sigma_Gamma,SSG,Multiple)
{
		log_p_Gamma_star=rep(0,Multiple)
		log_p_Gamma_2star=rep(0,Multiple)
		p=rep(0,Multiple)
		
		#if(det(sigma) == 0) {print(sigma); break}
		
		gamma_opt <-optim(Gamma,Z=Z,X=X,SG=Sigma_Gamma,phi=Phi,log_p_gamma,method="Nelder-Mead",
				hessian=TRUE,control=list(maxit=10,fnscale=-1))
		
		gamma=gamma_opt$par
		sigma=solve(-gamma_opt$hessian)
		diag(sigma)=diag(sigma)*SSG
		
		Gamma_star=t(rmvt(n=Multiple,gamma,sigma=sigma,df=df))
			
			for (i in 1:Multiple){
				log_p_Gamma_star[i]=log_p_gamma(Gamma_star[,i],Z,X,Sigma_Gamma,Phi)

			}
			
			
			#control overfloat, -max(log_p_Gamma_star);
			p=exp(log_p_Gamma_star-max(log_p_Gamma_star))/sum(exp(log_p_Gamma_star - max(log_p_Gamma_star)))
			
			#in case there is still overfloat;
			p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
			j=sample(c(1:Multiple),1,prob=p);	
							
			
			Gamma_2star=t(rmvt(n=Multiple-1,Gamma_star[,j],sigma=sigma,df=df))
			Gamma_2star <-cbind(Gamma_2star,Gamma)
			
			for (i in 1:Multiple){
				log_p_Gamma_2star[i]=log_p_gamma(Gamma_2star[,i],Z,X,Sigma_Gamma,Phi)
			}
			
											
			#control overfloat
			num=sum(exp(log_p_Gamma_star -max(log_p_Gamma_star)))
			den=sum(exp(log_p_Gamma_2star -max(log_p_Gamma_star)))
											
			rho=min(1,num/den)
			
			#in case overfloat again
			if(is.na(rho)) {rho=0.5};
			
			
			accp=0;
			u=runif(1)
			if(u<rho){
				Gamma=Gamma_star[,j]
				accp=1;
			
			}
		
		
		return(list(par=Gamma, accp=accp))
}

