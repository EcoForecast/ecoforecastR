#' GLOBAL UNCERTAINTY PARTITIONING
#' 
#' Function to calculate Sobol indices (relative variances) from 12 hindcast experiments
#' For each experiment, provide a matrix where rows = ensembles and columns = dates
#' 
#' @param A All B, as in: forecastN(ICB,XB,thetaB,alphaB,QmcB,zpB,Nmc) 
#' @param B All A, as in: forecastN(IC,X,theta,alpha,Qmc,zp,Nmc) 
#' @param ABp A with B’s parameter samples
#' @param ABi A with B’s initial condition samples
#' et cetera
#' output: list of two matrices (relative and total variances)
#'  
calculate_sobol_indices <- function(A, B, ABp, ABi, ABd,
																		ABz, ABpi, ABpd, ABpz,
																		ABid, ABiz, ABdz){
	
	Nf = B
	Nf.ci = apply(Nf,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
	
	## Saltelli eq 4.6: f0 = E(Y)
	mu = matrix(apply(Nf,2,mean,na.rm=TRUE),nrow(Nf),ncol(Nf),byrow = TRUE)
	
	# Subtract mu from each forecast experiment
	## Saltelli eq 4.7: fi = E(Y|Xi)-E(Y)
	A = A - mu
	B = B - mu
	
	# Main effects
	AB= list()
	AB[[1]] = ABi - mu
	AB[[2]] = ABd - mu
	AB[[3]] = ABp - mu
	AB[[4]] = ABz - mu
	names(AB)=c("i","d","p","z")
	
	# Interaction terms
	AB[[5]] = ABid - mu
	AB[[6]] = ABpi - mu
	AB[[7]] = ABiz - mu
	AB[[8]] = ABpd - mu
	AB[[9]] = ABdz - mu
	AB[[10]]= ABpz - mu
	names(AB)[5:10]=c("id","ip","iz","dp","dz","pz")
	for(i in 5:10){
		AB[[i]] = AB[[i]] - AB[[substr(names(AB)[i],1,1)]] - AB[[substr(names(AB)[i],2,2)]]
	}
	
	
	NT = ncol(A)
	d = length(AB)
	f0f0 = apply(A*B,2,mean,na.rm=TRUE)
	f02 = apply(A,2,mean,na.rm=TRUE)^2 ## Saltelli eqn 4.22
	VE = EV = matrix(NA,d,NT)
	
	## main effects: Si = V(fi)/V(Y)## main effects: Si = V(fi)/V(Y)
	for(i in 1:d){
		VEi = A*AB[[i]]   ## version from Saltelli
		VE[i,] = apply(VEi,2,mean,na.rm=TRUE) - f0f0 ## Var(E(Y|Xi)); f0f0 from Rcpp
	}
	S = t(t(VE)/apply(VE,2, sum))
	
	## interactions: Sij = Vij/V(Y)
	## Saltelli 4.10c: Vij =V(fij(Xi,Xj))=V(E(Y|Xi,Xj))−V(E(Y|Xi))−V(E(Y|Xj))
	
	## Total effects
	for(i in 1:d){
		EVi = B*AB[[i]]
		EV[i,] = apply(EVi,2,mean,na.rm=TRUE) - f0f0 ## E(Var(Y|Xi))
	}
	ST = t(1-t(EV)/apply(VE,2,sum))
	
	rownames(S) <- names(AB)
	rownames(ST) <- names(AB)
	
	out <- list("Relative effects" = S, "Total effects" = ST)
	return(out)
}

##### EXAMPLE: run function on each cyanobacteria forecast ------

#devtools::install_github("eco4cast/EFIstandards")
#devtools::install_github("EcoForecast/ecoforecastR")
#install.packages("jqr", configure.vars=c("LIB_DIR='/share/pkg.7/libjq/1.6/install/lib' INCLUDE_DIR='/share/pkg.7/libjq/1.6/install/include'"))
pacman::p_load(ncdf4, tidyverse, EML, emld, uuid, EFIstandards, ecoforecastR, stringr)

# Read in necessary helper functions
source('https://raw.githubusercontent.com/melofton/Bayes_forecast_WG/NEFI_UA_synthesis/0_Function_library/metadata_functions.R') 

# Another helper function
pivot_to_matrix <- function(tibble) {
	out_matrix <- pivot_wider(tibble, 
														names_from = forecast_valid_time, 
														values_from = Gloeo_abundance) %>% 
		select(-c(forecast_issue_time, depth, lon, lat, ensemble)) %>% 
		as.matrix()
	return(out_matrix)
}


# Get date list from directory contents
file.list <- list.files(path = "/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/",
												pattern = "cyanobacteria")
date_list <- str_extract(file.list, '(?<=-).*(?=-)') %>% unique()
date_list = "20130523" # for testing



##### Loop through all dates -------

# For outputs
relative_effects_list <- list()
total_effects_list <- list()

for (issue_date in date_list){
	
	# Read in all forecasts with given issue date, and pivot to a matrix where rows = ensembles and columns = dates
	A_orig = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-A.nc"),
													 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	B_orig    = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-B.nc"),
															var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABp  = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABp.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABi  = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABi.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABd  = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABd.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABz  = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABz.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABpi = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABpi.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABpd = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABpd.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABpz = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABpz.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABid = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABid.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABiz = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABiz.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	ABdz = nc_hindcast_get(paste0("/projectnb/dietzelab/NEFI_SOBOL/cyanobacteria/cyanobacteria-", issue_date, "-ABdz.nc"),
												 var_name = 'Gloeo_abundance') %>% pivot_to_matrix()
	
	## Flip A and B to match Saltelli definition of Ci = all B except Ai
	B = A_orig
	A = B_orig
	
	ind_out <- calculate_sobol_indices(A, B, ABp, ABi, ABd,
																		 ABz, ABpi, ABpd, ABpz,
																		 ABid, ABiz, ABdz)
	# Save output to list
	relative_effects_list[[issue_date]] <- ind_out[['Relative effects']]
	total_effects_list[[issue_date]]    <- ind_out[['Total effects']]
	
}

##### VISUALIZE
date_S <- relative_effects_list[[1]]
NT = ncol(date_S)
N.cols <- c("black","green","red","orange","cyan","blue") ## set colors
N.cols = RColorBrewer::brewer.pal(n=nrow(date_S),name = "Set3") ## main effects
date_Sc = apply(date_S,2,cumsum)
plot(1:NT,date_Sc[1,],ylim=c(0,1),type='n',main="Relative Variance: In-Sample",ylab="Proportion of Variance",xlab="time")
ciEnvelope(1:NT,rep(0,ncol(date_Sc)),date_Sc[1,],col=N.cols[1])
for(i in 2:nrow(date_Sc)){
	ciEnvelope(1:NT,date_Sc[i-1,],date_Sc[i,],col=N.cols[i])
}
legend("topleft",legend = rev(rownames(date_Sc)),col=rev(N.cols),lty=1,lwd=5) 
