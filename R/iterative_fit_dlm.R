##' @name iterative_fit_dlm
##' @title iterative_fit_dlm
##' @author Mike Dietze
##' @export
##' @param model list containing the following elements
##' \itemize{
##'  \item{obs}{ column name of the observed data. REQUIRED}
##'  \item{fixed}{ formula for fixed effects. Response variable is optional but should be 'x' if included}
##'  \item{random}{ not implemented yet; will be formula for random effects}
##'  \item{n.iter}{ number of mcmc iterations}
##' }
##' @param data  data frame containing observations and covariates
##' @param method  "increment" to update an analysis with a new increment of data (i.e. Analysis step of data assimilation); "refit" to fit from t=0; or "iterative" to fit just next increment using previous posterior as prior
##' @param dt    number of time steps fit in each iteration, ignored for "increment"
##' @param nf    number of time steps to forecast into the future
##' @param prev  output from previous iteration (only used when method=increment; `model` priors will be ignored)
##' @description Iteratively fits a Bayesian state-space dynamic linear model using JAGS
iterative_fit_dlm <- function(model = NULL,
                              data,
                              method = "increment",
                              dt = 10,
                              nf = 1,
                              prev = NULL) 
{
  if(FALSE){ ## settings for testing
    model=list(obs="Y",fixed=NULL)
    data=data.frame(Y=Y)
    method = "refit"
    dt = 10
    nf = 1
  }
  
  ## method
  method = tolower(method)
  if(!(method %in% c("refit"))){
    print(paste("Method =",method,"not supported yet"))
    return(NULL)
  }
  
  ## data
  data  <- as.data.frame(data)
  nstep <- nrow(data) / dt
  
  ## state observations
  obs   <- model$obs
  if(is.null(obs)){
    print("Observations not included in model. Please add the variable 'obs' to the model list")
  } else {
    if(obs %in% names(data)){
      OBS = data[,obs]
    } else {
      print(paste("Could not find",obs,"in the provided data frame"))
      return(NULL)
    }
  }
    
  ## set up priors
  tau_a <- tau_r <- matrix(0.1, nstep + 1,2)
  fixed = model$fixed
  if(fixed == ""){
    nbeta = 0
  } else {
   # nbeta <- ncol(BuildZ(fixed,data))
    fixedX <- sub("~","",fixed, fixed=TRUE)
    lm.terms <- gsub("[[:space:]]", "", strsplit(fixedX,split = "+",fixed=TRUE)[[1]])  ## split on + and remove whitespace
    lm.terms <- unlist(strsplit(lm.terms,split = "-",fixed=TRUE))  ## split on -
    nbeta <- length(lm.terms)
  }
  if(!is.null(nbeta) & nbeta > 0){
    betaCOV <- array(NA,c(nstep+1,nbeta,nbeta))
    for(i in seq_len(nstep)) betaCOV[i,,] <- solve(diag(0.001,nbeta))
  }
  
  ## set up storage
  paramStats <- list()
  xp <- array(NA, c(nstep,6,nf)) ## prediction lowerCI, median, upperCI, mean, var, observation quantile
  dicI <- matrix(NA, nstep, 2)   ## Iterative updating of DIC
  for (i in seq_len(nstep)) {
      ## prep data, pad with NA for prediction
      if(method %in% c("refit","increment")){
        start = 1
      } else {
        start = 1 + (i - 1) * dt
      }
      stop  = i*dt
      MyDAT <- data[start:(stop+nf),,drop=FALSE]
      if(nf > 0){
        MyDAT[(stop+1:nf),obs] <- NA ## forecast by setting obs to missing data
      }

      ## fit model
      mcI <- fit_dlm(model=model,data=MyDAT)

      ## calculate DIC
      dicI[i, 1] <- sum(mcI$DIC$deviance)
      dicI[i, 2] <- sum(mcI$DIC$penalty)
      
      ## summerize parameters
      paramStats[[i]] <- summary(mcI$params)
      
      ## update BETA priors
      mcI.m  <- as.matrix(mcI$params) ## convert to matrix, remove burn-in
      betaID <- which(grepl(pattern = "^beta",colnames(mcI.m)))
      if(length(betaID)>0){
        if(length(betaID)>1){
          betaCOV[i+1,,] <- cov(mcI.m[,betaID])
        } else {
          betaCOV[i+1,,] <- var(mcI.m[,betaID])
        }
      } else { 
        betaCOV <- NULL  ## nbeta = 0
      }

      ## update TAU priors
      tauID  <- which(grepl(pattern = "^tau",colnames(mcI.m)))
      m <- paramStats[[i]]$statistics[tauID,1]
      s <- paramStats[[i]]$statistics[tauID,2]
      tau_a[i + 1,] <- m ^ 2 / s^2
      tau_r[i + 1,] <- m / s^2
      
      ## save predictions
      if (nf > 0) {
        predID <- stop + 1:nf
        mcI.x  <- as.matrix(mcI$predict[,predID])
        xp[i,1:3,] <- apply(mcI.x,2,quantile,c(0.025, 0.5, 0.975)) ## prediction quantiles
        xp[i,4,]   <- apply(mcI.x,2,mean)                          ## prediction mean
        xp[i,5,]   <- apply(mcI.x,2,var)                           ## prediction variance
        for(j in seq_along(predID)){
          xp[i,6,j]  <- sum(mcI.x[,j] < OBS[stop + j]) / nrow(mcI.m) ## Bayes pval, predictive error
        }
      }
  } ## end loop over steps
  
  return(list(paramStats=paramStats,xp=xp,dic=dicI,
              tau_a=tau_a,tau_r=tau_r,betaCOV=betaCOV,
              method=method,dt=dt,nf=nf))
}
