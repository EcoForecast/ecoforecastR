##' @name fit_dlm
##' @title fit_dlm
##' @author Mike Dietze
##' @export
##' @param model list containing the following elements
##' \itemize{
##'  \item{obs}{column name of the observed data. REQUIRED}
##'  \item{fixed}{formula for fixed effects. Response variable is optional but should be 'x' if included}
##'  \item{random}{not implemented yet; will be formula for random effects}
##'  \item{n.iter}{number of mcmc iterations}
##' }
##' @param data  data frame containing observations and covariates
##' @param dic   whether or not to calculate DIC
##' @description Fits a Bayesian state-space dynamic linear model using JAGS
fit_dlm <- function(model=NULL,data,dic=TRUE){

  obs    = model$obs
  fixed  = model$fixed
  random = model$random
  n.iter = ifelse(is.null(model$n.iter),5000,model$n.iter)

  data = as.data.frame(data)
  
  out.variables = c("x","tau_obs","tau_add")
  Pformula = NULL
  
  ## observation design matrix
  if(is.null(obs)){
    print("Observations not included in model. Please add the variable 'obs' to the model list")
  } else {
    if(length(grep("~",obs)) == 0){ ## no formula, assuming obs is just a variable
      if(obs %in% names(data)){
        OBS = data[,obs]
      } else {
        print(paste("Could not find",obs,"in the provided data frame"))
        return(NULL)
      }
    } else {  ## obs is a formula
      print("obs formulas not implemented yet")
      return(NULL)
    }
  }
  
  #### prep data
  mydat<-list(OBS=OBS,n=length(OBS),x_ic = 0,tau_ic = 0.00001,a_obs=0.1,r_obs=0.1,a_add=0.1,r_add=0.1)
  
  
  ## process design matrix
  if(is.null(fixed) | fixed == ""){
    fixed = NULL
  } else {
    if(is.null(data)) print("formula provided but covariate data is absent:",fixed)
    design <- ParseFixed(fixed,data,
                         update=list(out.variables=out.variables,
                                     data = mydat))
#    Z = as.matrix(Z[,-which(colnames(Z)=="(Intercept)")])
    if(sum(is.na(design$data))>0){
      print("WARNING: missing covariate data")
      print(apply(is.na(design$data),2,sum))
    }
  }
  ## alternatively might be able to get fixed and random effects simultaneously using
  ## lme4::lFormula(formula("x ~ FIXED + (1|FACTOR)"),na.action=na.pass)
  ## e.g. foo = lme4::lFormula(formula("x ~ PAR + (1+PAR|DOY)"),na.action = na.pass)

  
  
  #### Define JAGS model
  my.model = "  
  model{
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)

  #### Random Effects
  #RANDOM  tau_alpha~dgamma(0.1,0.1)
  #RANDOM  for(i in 1:nrep){                  
  #RANDOM   alpha[i]~dnorm(0,tau_alpha)
  #RANDOM  }

  #### Fixed Effects
  ##BETAs
  ##MISSING_MU
  
  #### Data Model
  for(t in 1:n){
    OBS[t] ~ dnorm(x[t],tau_obs)
    ##MISSING
  }
  
  #### Process Model
  for(t in 2:n){
    mu[t] <- x[t-1] ##PROCESS
    x[t]~dnorm(mu[t],tau_add)
  }

  }"
  
 
  #### prep model
  if(!is.null(fixed)){
    
    ## Insert regression priors
    my.model = sub(pattern="##BETAs",design$Xpriors,my.model)  
    out.variables = design$out.variables
    mydat = design$data
    Pformula = design$Pformula
    Pnames = unique(design$Pnames)

    ## missing data model
    if(!is.null(design$MDprior)){
      my.model <- sub(pattern="##MISSING_MU",design$MDprior,my.model)
      my.model <- sub(pattern="##MISSING",design$MDformula,my.model)
    }
  }
  
  ## RANDOM EFFECTS
  if(!is.null(random)){
    my.model = gsub(pattern="#RANDOM"," ",my.model)
    out.variables = c(out.variables,"tau_alpha","alpha")  
    Pformula = " + alpha[rep[i]]"
    ## *****
    ## need to do work here to specify indicator variables for random effects explictly
    ## *****
  }
  
  if(!is.null(Pformula)) my.model = sub(pattern="##PROCESS",Pformula,my.model)
  
  ## Define initial conditions

  
  print(my.model)
  ## initialize model
  mc3 <- rjags::jags.model(file=textConnection(my.model),data=mydat,
                    n.chains=3)
  
  mc3.out <- rjags::coda.samples(model=mc3, variable.names=out.variables, n.iter=n.iter)
  
  ## split output
  out = list(params=NULL,predict=NULL,model=my.model,data=mydat)
  mfit = as.matrix(mc3.out,chains=TRUE)
  pred.cols = union(grep("x[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
  chain.col = which(colnames(mfit)=="CHAIN")
  out$predict = mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
  out$params   = mat2mcmc.list(mfit[,-pred.cols])
  if(dic) out$DIC <- dic.samples(mc3, 1000)
  return(out)
  
}  ## end fit_dlm



BuildZ <- function(fixed, data) {
  if (toupper(fixed) == "RW") {
    return(NULL)
  } else {
    fixed = ifelse(length(grep("~", fixed)) == 0, paste("~", fixed), fixed)
    fixed = sub("x*~", "~", x = fixed)
    options(na.action = na.pass)
    #  Z = with(data,model.matrix(formula(fixed),na.action=na.pass))
    Z = model.matrix(formula(fixed),
                     data = model.frame(formula(fixed), data),
                     na.action = na.pass)
    return(Z)
  }
}