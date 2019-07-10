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
  a_add <- r_add <- rep(0.1, nstep + 1)
  
  ## set up storage
  mp <- xp <- tauI <- matrix(NA, nstep, 3) ## median and CI
  mq <- xq <- rep(NA, nstep) #quantile of prediction
  mv <- xv <- rep(NA, nstep) #variance of prediction
  dicI <- matrix(NA, nstep, 2)     #Iterative updating of DIC
    
  for (i in seq_len(nstep)) {
      ## prep data, pad with NA for prediction
      if(method %in% c("refit","increment")){
        start = 1
      } else {
        start = 1 + (i - 1) * dt
      }
      stop  = i*dt
      mydat <- data[start:(stop+nf),,drop=FALSE]
      if(nf > 0){
        mydat[(stop+1:nf),obs] <- NA ## forecast by setting obs to missing data
      }

      ## fit model
      mcI <- fit_dlm(model=model,data=data)

      ## calculate DIC
      DIC <- dic.samples(mcI, 1000)
      dicI[i, 1] <- sum(DIC$deviance)
      dicI[i, 2] <- sum(DIC$penalty)
      
      ## save parameters
      mcI.m <- as.matrix(mcI.out) ## convert to matrix, remove burn-in
      tauI[i, ] <- quantile(1 / sqrt(mcI.m[, 1]), c(0.025, 0.5, 0.975))
      tbar <- mean(1 / sqrt(mcI.m[, 1]))
      
      ## update priors
      m  = mean(mcI.m[, 1])
      s2 = var(mcI.m[, 1])
      a_add[i + 1] <- m ^ 2 / s2
      r_add[i + 1] <- m / s2
      
      if (FALSE) {
        hist(mcI.m, probability = TRUE)
        xseq <- seq(min(mcI.m), max(mcI.m), length = 1000)
        lines(xseq, dgamma(xseq, a_add[i + 1], r_add[i + 1]))
      }
      
      ## save predictions
      xp[i, ] <-
        quantile(mcI.m[, 2], c(0.025, 0.5, 0.975)) ## prediction quantiles
      xq[i]  <-
        sum(mcI.m[, 2] < Y[i * dt + 1]) / nrow(mcI.m) ## Bayes pval, predictive error
      xv[i]  <-
        var(mcI.m[, 2])                         ## prediction variance
      
      ## prediction @ parameter mean
      mp[i, ] <- Y[i * dt] + c(-1.96, 0, 1.96) * tbar
      mq[i]  <- pnorm(Y[i * dt + 1], Y[i * dt], tbar)
      mv[i]  <- tbar ^ 2
    }
}
