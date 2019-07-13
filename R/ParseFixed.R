ParseFixed <- function(fixed,cov.data,update=NULL,ancillary.dims=NULL){
  
  if(FALSE){
    ## DEV TESTING FOR X, polynomial X, and X interactions
    fixed <- "X + X^3 + X*bob + bob + dia + X*Tmin[t]" ## faux model, just for testing jags code
  }
  
  ## set up string variables (OK for many to start NULL)
  data = update$data
  out.variables = update$out.variables
  Pformula = update$Pformula
  Xpriors = update$Xpriors
  Pnames = update$Pnames
  MDprior = update$MDprior
  MDformula = update$MDformula
  
  ## parse if working with a single time series or additional dimensions
  if(is.null(ancillary.dims)){
    AD=""
  } else {
    AD=ancillary.dims
  }
  
  ## Design matrix
  if (is.null(fixed)) {
    Xf <- NULL
  } else {
    
    ## check for covariate data (note: will falsely fail if only effect is X)
    if (is.null(cov.data)) {
      print("formula provided but covariate data is absent:", fixed)
    } else {
      cov.data <- as.data.frame(cov.data)
    }
    
    ## check if there's a tilda in the formula
    if (length(grep("~", fixed)) == 0) {
      fixed <- paste("~", fixed)
    }
    
    ## First deal with endogenous terms (X and X*cov interactions)
    fixedX <- gsub("[[:space:]]", "", sub("~","",fixed, fixed=TRUE))
    lm.terms <- unlist(strsplit(fixedX,split = "-",fixed=TRUE))  ## split on -
    if(lm.terms[1] == ""){ ## was negative
      lm.terms = lm.terms[-1]
      lm.terms[1] = paste0("-",lm.terms[1])
     }
    if(length(lm.terms)>1){lm.terms[2] = paste0("-",lm.terms[2])} ## restore later minus
    lm.terms <- unlist(strsplit(lm.terms,split = "+",fixed=TRUE))  ## split on + and remove whitespace
    X.terms <- strsplit(lm.terms,split = c("^"),fixed = TRUE)
    X.terms <- sapply(X.terms,function(str){unlist(strsplit(str,,split="*",fixed=TRUE))})
    X.terms <- which(sapply(X.terms,function(x){any(toupper(x) == "X")}))
    if(length(X.terms) > 0){
      ## rebuild fixed without X.terms
      fixed <- paste("~",paste(lm.terms[-X.terms],collapse = " + "))  
      
      ## isolate terms with X
      X.terms <- lm.terms[X.terms]
      for(i in seq_along(X.terms)){
        
        myBeta <- NULL
        Xformula <- NULL
        if(length(grep("*",X.terms[i],fixed = TRUE)) == 1){  ## INTERACTION
          
          myIndex <- "[t-1]"             ### changed this from i to t, may break things 7/10/19 ***
          covX <- strsplit(X.terms[i],"*",fixed=TRUE)[[1]] 
          covX <- covX[-which(toupper(covX)=="X")] ## remove X from terms
          
          ##is covariate fixed or time varying?
          tvar <-  length(grep("[t]",covX,fixed=TRUE)) > 0           
          if(tvar){
            covX <- sub("[t]","",covX,fixed = TRUE)
            if(!(covX %in% names(data))){
              ## add cov variables to data object
              data[[covX]] <- time_data[[covX]]
            }
            check.dup.data(data,"covX")
            
            myIndex <- "[i,t]"
          } else {
            ## variable is fixed
            if(covX %in% colnames(cov.data)){ ## covariate present
              if(!(covX %in% names(data))){
                ## add cov variables to data object
                data[[covX]] <- cov.data[,covX]
              }
              check.dup.data(data,"covX2")
              
            } else {
              ## covariate absent
              warning("covariate absent from covariate data:", covX)
            }
            
          } ## end fixed or time varying
          
          myBeta <- paste0("betaX_",covX)
          Pnames = c(Pnames,covX)
          Xformula <- paste0(myBeta,"*x[",AD,"t-1]*",covX,myIndex)  ## was x[i,t-1]
          
        } else if(length(grep("^",X.terms[i],fixed=TRUE))==1){  ## POLYNOMIAL
          powX <- strsplit(X.terms[i],"^",fixed=TRUE)[[1]] 
          powX <- powX[-which(toupper(powX)=="X")] ## remove X from terms
          myBeta <- paste0("betaX",powX)
          Xformula <- paste0(myBeta,"*x[",AD,"t-1]^",powX)
          
        } else {  ## JUST X
          myBeta <- "betaX"
          Xformula <- paste0(myBeta,"*x[",AD,"t-1]")
        }
        
        Pformula <- paste(Pformula,"+",Xformula)
        
        ## add priors
        Xpriors <- paste(Xpriors,"     ",myBeta,"~dnorm(0,0.001)\n")
        
        ## add to out.variables
        out.variables <- c(out.variables, myBeta)
        
      }  ## END LOOP OVER X TERMS
      
    }  ## end processing of X terms
    Pnames = unique(Pnames)
    
    ############ build DESIGN MATRIX from formula #########
    fixedX <- sub("~","",fixed, fixed=TRUE)
    lm.terms <- gsub("[[:space:]]", "", strsplit(fixedX,split = "+",fixed=TRUE)[[1]])  ## split on
    if(lm.terms[1] != ""){
      Xf = model.matrix(formula(fixed),
                        data = model.frame(formula(fixed), cov.data,na.action = na.pass),
                        na.action = na.pass)
      Xf.cols <- colnames(Xf)
      Xf.cols <- sub(":","_",Xf.cols) ## for interaction terms, switch separator
      colnames(Xf) <- Xf.cols
      # Xf.cols <- Xf.cols[Xf.cols != "(Intercept)"]
      # Xf      <- as.matrix(Xf[, Xf.cols])
      # colnames(Xf) <- Xf.cols
      ##Center the covariate data
      #    Xf.center <- apply(Xf, 2, mean, na.rm = TRUE)
      #    Xf      <- t(t(Xf) - Xf.center)
    } else {Xf <- NULL} ## end fixed effects parsing
    
    ## drop -1 term, isn't part of design so shouldn't get a beta
    if(ncol(Xf) == 0) Xf <- NULL
    
    ## build formula in JAGS syntax
    if (!is.null(Xf)) {
      Xf.names <- gsub(" ", "_", colnames(Xf))  ## JAGS doesn't like spaces in variable names
      Xf.names <- gsub("(", "", Xf.names,fixed=TRUE)  ## JAGS doesn't like parentheses in variable names
      Xf.names <- gsub(")", "", Xf.names,fixed=TRUE)
      
      ## remove items from design matrix that are already in model
      real.names <- Xf.names[Xf.names %in% Pnames]
      sel <- which(Xf.names %in% real.names)
      if(length(sel)>0){
        Xf <- as.data.frame(Xf)
        Xf <- Xf[,-sel,drop=FALSE]
        Xf.names <- Xf.names[-sel]
      }
      
      ## append to process model formula: Xf
      if(ncol(Xf)>0){
        Pformula <- paste(Pformula,
                        paste0("+ beta", Xf.names, "*Xf[t,", seq_along(Xf.names), "]", collapse = " "))  # was Xf[rep[i]
        Xpriors <- paste(Xpriors,paste0("     beta", Xf.names, "~dnorm(0,0.001)", collapse = "\n"))
        MDprior <- paste(MDprior,
                         "for(j in 1:",ncol(Xf),"){\n",
                         "   muXf[j] ~ dnorm(0,0.001)\n",
                         "   tauXf[j] ~ dgamma(0.01,0.01)\n",
                         "}\n")
        MDformula <- paste(MDformula,
                           paste0("Xf[t,",seq_along(Xf.names),
                                  "] ~ dnorm(muXf[",seq_along(Xf.names),
                                  "],tauXf[",seq_along(Xf.names),"])",collapse="\n")
                          )
        out.variables <- c(out.variables, paste0("beta", Xf.names))
      }
      ## append using real names
      if(length(real.names) > 0){
        Pformula <- paste(Pformula,
                          paste0("+ beta", real.names, "*",real.names,"[t]", collapse = " "))
        Xpriors <- paste(Xpriors,paste0("     beta", real.names , " ~ dnorm(0,0.001)", collapse = "\n"))
        out.variables <- c(out.variables, paste0("beta", real.names))
      }
      
      ## create 'rep' variable if not defined
#     if(is.null(data$rep)){
#        data$rep <- seq_len(nrow(Xf))
#      }

            ## update variables for JAGS to track
      data[["Xf"]] <- Xf
    }
    ## missing data model for Pnames (do only once across both interactions and Xf)
    missCol <- which(Pnames != "Intercept")
    if(length(missCol)>0){
      Pmiss <- Pnames[missCol]
      MDprior <- paste(
        paste0("mu",Pmiss,"~dnorm(0,0.001)",collapse="\n"),"\n",
        paste0(" tau",Pmiss,"~dgamma(0.01,0.01)",collapse="\n")
      )
      MDformula <- paste0(Pmiss,"[t] ~ dnorm(mu",Pmiss,",tau",Pmiss,")",collapse="\n")
    }
    
    check.dup.data(data,"Xf")
  
  } ## END FIXED IS NOT NULL
  
  return(list(Pformula=Pformula,out.variables=out.variables,Xpriors=Xpriors,
              MDprior=MDprior, MDformula=MDformula,data=data))
}

check.dup.data <- function(data,loc){
  if(any(duplicated(names(data)))){warning("duplicated variable at ",loc," ",names(data))}
}

if(FALSE){
  ## DUMPING GROUND
  TreeDataFusionMV <- sub(pattern = "## ENDOGENOUS BETAS", Xpriors, TreeDataFusionMV)
  TreeDataFusionMV <- sub(pattern = "## FIXED EFFECTS BETAS", Xf.priors, TreeDataFusionMV)
  
}