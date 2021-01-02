cmlFDR_GammaDist=function (Z,X,nIter=1100,burnIn=100,thin=5,initNULL=0.95,simulate=FALSE,
			SSA=1,SSG=1,MA=3,MG=3,mu=0.68,theoNULL=FALSE,inits=NULL)
{
#SSA:scale the diagnal of the covariance matrix in MH of Alpha draw to increase/decrease the step size
#SSG:scale the diagnal of the covariance matrix in MH of Gamma draw to increase/decrease the step size


###########################
#Load packages/functions
###########################
	
	library(tmvtnorm)
	library(mnormt)
	library(magic)
	library(locfdr)	
	library(arm)
	library(pscl)

	source("fcns.R")

	if(colnames(X)[1] != "Intcpt") X=cbind(Intcpt=1,X)
	#print(X[1:10,])
	# hyperparameters
	
	B01=0.001;B02=0.001; #Prior: P(Beta) ~ Gamma(B01,B02)
	a0=0.001;b0=0.001; #Prior: P(sigma^2) ~ IVG(a0,b0) #b0 is scale;
	Sigma_Gamma=matrix(0,nrow=dim(X)[2],ncol=dim(X)[2])
	diag(Sigma_Gamma)=10000; #Prior: P(Gamma) ~ N(0,Sigma_Gamma)

	Sigma_Alpha=matrix(0,nrow=dim(X)[2],ncol=dim(X)[2])
	diag(Sigma_Alpha)=10000; #Prior: P(Alpha) ~ N(0,Sigma_Alpha)
	df=4 # degrees of freedom for multivariate t proposal
	
	# data
	N=dim(X)[1]
	M=dim(X)[2]

	# Parameter arrays	
	
	ALPHA_array=list()
	BETA_array=list()
	if(theoNULL==FALSE){
		SIGMA_SQ_array=list()
	}
	GAMMA_array=list()
	Accp_Rate_array=list()  	#Alpha draw accept rate
	Accp_Rate_array_g=list()	#Gamma draw accept rate
	
	array_ind=0

	# Initialize parameters 

	if(is.null(inits)){
		Alpha=cbind(rep(0,M));Alpha_mean=Alpha
		pi0=initNULL #min(.98,locfdr(Z,nulltype=1,plot=0)$fp0[5,3]);
		gamma0=log((1-pi0)/pi0) #intercept for Non-NULL gamma
		Gamma=array(0,dim=c(M,1));Gamma_mean=Gamma; Gamma[1,]=gamma0
		Beta=0.1;Beta_mean=Beta
		Phi=1-as.numeric(abs(Z)< sort(abs(Z))[round(pi0*N)]);Phi_mean=0*Phi
		Sigma_sq=1
	}
	if(!is.null(inits)){
		last=length(inits[[1]])
		Alpha=cbind(inits[[1]][[last]]);Alpha_mean=Alpha
		Beta=inits[[2]][[last]];Beta_mean=Beta
		Gamma=cbind(inits[[3]][[last]]);Gamma_mean=Gamma
		if(theoNULL==FALSE){
			Sigma_sq=inits[[4]][[last]];Sigma_sq_mean=Sigma_sq
		}
		if(theoNULL==TRUE){
			Sigma_sq=1
		}
		Phi=inits[[5]];Phi_mean=0*Phi
	}		
	PHI_match_rate=NULL


	for(iter in 1:nIter){
		
#			print(iter)
		
			Z1=abs(Z[!Phi==0])
			X1=X[!Phi==0,]

	
		## Draw ALPHA 
			
			if(det(t(X1)%*%X1) != 0){ #avoid singular	
				obj=Draw_Alpha_log_M_mu(Alpha,Z1,X1,Beta,Phi,df,MA,Sigma_Alpha,SSA,mu)
				Alpha=obj$par;
			}
#			print("Alpha");
#			print(Alpha);
			
			
		## Draw BETA
		
			Beta=rgamma(1,shape=B01+sum(exp(X1%*%Alpha)),rate=B02+sum(Z1-mu))
#			print(paste("Beta",Beta));
			
					
		## Draw GAMMA
			
			#Gamma draw: Multiple try MH
			objg=Draw_Gamma_log_M(Gamma,Z,X,Phi,df,Sigma_Gamma,SSG,MG)

			Gamma=objg$par;
#			print("Gamma");
#			print(Gamma);
		

		if(theoNULL==FALSE){
		## Draw SIGMA_SQ;
		   	Z0=abs(Z[Phi==0])
			Sigma_sq=rigamma(1,alpha=a0+length(Z0)/2,beta=b0+(Z0%*%Z0)/2);
#			print(paste("Sigma_sq",Sigma_sq));
		}


		## Draw PHI; 	
		
			log_P_phi=cbind(X%*%Gamma,0)
			log_P_phi[abs(Z)>mu,1]=log_P_phi[abs(Z)>mu,1]+((exp(X[abs(Z)>mu,]%*%Alpha)-1)*
				log(abs(Z[abs(Z)>mu])-mu)-Beta*(abs(Z[abs(Z)>mu])-mu))+(exp(X[abs(Z)>mu,]%*%Alpha))*
				log(Beta)-lgamma(exp(X[abs(Z)>mu,]%*%Alpha))
			#log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi)-0.5*Z^2;
			
			if(theoNULL==TRUE){
				log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi)-0.5*Z^2;
			}else{
				log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi*Sigma_sq)-(1/(2*Sigma_sq))*Z^2;
			}

			P_phi=exp(log_P_phi)/apply(exp(log_P_phi),1,sum) #declare the variable P_phi;
			P_phi[abs(Z)>mu,1]=1/(1+exp(log_P_phi[abs(Z)>mu,2]-log_P_phi[abs(Z)>mu,1]))
			P_phi[,2]=1/(1+exp(log_P_phi[,1]-log_P_phi[,2]))

			P_phi[abs(Z)<=mu,1]=0;P_phi[abs(Z)<=mu,2]=1;
			
			Phi_new=Phi #declare variable Phi_new, create a vector;
			for(i in 1:N){
				Phi_new[i]=sample(c(1,0),size=1,replace=TRUE,prob=P_phi[i,])
			}
			Phi=Phi_new
			
#			if(simulate == TRUE) print(sum(Phi == Phi_true)/N)
			

		## Save results after thin
				
			if(iter%%thin==0 & iter>=burnIn){
				array_ind=array_ind+1
				ALPHA_array[[array_ind]]=Alpha
				
				GAMMA_array[[array_ind]]=Gamma

				BETA_array[[array_ind]]=Beta
				
				if(theoNULL==FALSE){
					SIGMA_SQ_array[[array_ind]]=Sigma_sq
				}

				if(det(t(X1)%*%X1) != 0){
					Accp_Rate_array[[array_ind]]=obj$accp
				}
				else Accp_Rate_array[[array_ind]]=0
				
				Accp_Rate_array_g[[array_ind]]=objg$accp;

#				print("Alpha mean:");
				Alpha_mean=((array_ind-1)*Alpha_mean+ALPHA_array[[array_ind]])/array_ind
#				if(simulate==TRUE) print(cbind(Alpha_mean,Alpha_true))
#				if(simulate==FALSE) print(Alpha_mean)
				
				
#				print("Beta mean:");
#				if(simulate==TRUE) print(cbind(mean(as.numeric(BETA_array)),BETA_true))
#				if(simulate==FALSE) print(mean(as.numeric(BETA_array)))
				
#				if(theoNULL==FALSE){
#					print("Sigma_sq mean:");
#					if(simulate==TRUE) print(cbind(mean(as.numeric(SIGMA_SQ_array)),SIGMA_SQ_true))
#					if(simulate==FALSE) print(mean(as.numeric(SIGMA_SQ_array)))

#				}

#				print("Gamma mean:");
				Gamma_mean=((array_ind-1)*Gamma_mean+GAMMA_array[[array_ind]])/array_ind
#				if(simulate==TRUE) print(cbind(Gamma_mean,Gamma_true))
#				if(simulate==FALSE) print(Gamma_mean)
				
#				print(paste("Multiple-try MH Accept Rate for Alpha (mean):",mean(as.numeric(Accp_Rate_array))));
#				print(paste("Multiple-try MH Accept Rate for Gamma (mean):",mean(as.numeric(Accp_Rate_array_g))));				
			
				if(simulate==TRUE) {
					PHI_match_rate=rbind(PHI_match_rate,sum(Phi==Phi_true)/N);
#					print(paste("Phi matching rate (mean):",mean(PHI_match_rate)))
				}
				
				#probability of each SNP being Non-NULL, average of Phi over the iterations saved;
				Phi_mean=((array_ind-1)*Phi_mean+Phi)/array_ind;
					
			
				results=list()
				results[[1]]=ALPHA_array
				results[[2]]=BETA_array
				results[[3]]=GAMMA_array
				if(theoNULL==TRUE){	
					results[[4]]=1
				}
				if(theoNULL==FALSE){
					results[[4]]=SIGMA_SQ_array
				}
				results[[5]]=Phi
				#save(X,Z,results,file="mcmc_intermediate_outputs.R")

			}

		

	}	

	#Calculate SD;
	#Alpah
#	print("Alpha SD:");
#	for(i in 1:dim(X)[2]){
#    		print(sd(as.numeric(unlist(lapply(ALPHA_array, function(x) x[i])))))
#	}

	#beta
#	print(paste("Beta SD:",sd(as.numeric(BETA_array))))

#	if(theoNULL==FALSE){
		#Sigma_sq
#		print(paste("Sigma SD:",sd(as.numeric(SIGMA_SQ_array))))
#	}

	#Gamma
#	print("Gamma SD:");
#	for(i in 1:dim(X)[2]){
#    		print(sd(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i])))))
#	}


	#return results;
	#return results;
	results=list()
	results[[1]]=ALPHA_array
	results[[2]]=BETA_array
	results[[3]]=GAMMA_array
	if(theoNULL==TRUE){	
		results[[4]]=1
	}
	if(theoNULL==FALSE){
		results[[4]]=SIGMA_SQ_array
	}
	results[[5]]=Phi


	return(results)
}
