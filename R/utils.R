##' @name mat2mcmc.list
##' @title mat2mcmc.list
##' @export
##' @author Mike Dietze
##' @description convert a matrix to a CODA mcmc.list
##' @param w  matrix
mat2mcmc.list <- function(w) {
  temp <- list()
  chain.col <- which(colnames(w) == "CHAIN")
  for (i in unique(w[, "CHAIN"])) {
    temp[[i]] <- coda:::as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
  }
  return(as.mcmc.list(temp))
}

##' @name ciEnvelope
##' @title ciEnvelope
##' @export 
##' @author Mike Dietze
##' @description plot a confidence/credible interval
##' @param x x-axis data
##' @param ylo y-axis lower confidence bound
##' @param yhi y-axis upper confidence bound
##' @param ... optional graphical parameters
ciEnvelope <- function(x, ylo, yhi, ...) {
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi), ylo[1])), 
    border = NA, ...)
}

##' add transparency to a color
##' @author Mike Dietze
##' @export
##' @param col color (by name, number, or hex string)
##' @param alpha transparancy on a 0-1 scale (1=opaque)
col.alpha <- function(col,alpha=1){
  rgb = col2rgb(col)
  rgb(rgb[1],rgb[2],rgb[3],alpha*255,maxColorValue=255)
}

##' weighted quantile
##' @author Mike Dietze
##' @export
##' @param x numeric vector whose quantiles are wanted
##' @param wt vector of weights
##' @param q  quantile (at the moment must be scalar) 
##' generalization of the base quantile function to allow values to be weighted. Useful for constructing interval estimates when using a particle filter or similar algorithm
wtd.quantile <- function(x,wt,q){ 
  ord <- order(x)
  wstar <- cumsum(wt[ord])/sum(wt)
  qi <- findInterval(q,wstar); qi[qi<1]=1;qi[qi>length(x)]=length(x)
  return(x[ord[qi]])
}

##' @name plot_ss
##' @title plot_ss
##' @export
##' @author Mike Dietze
##' @description plot time series output of a state-space model
##' @param time time axis values
##' @param fit fit_dlm output object
##' @param add add lines to current plot. Default is FALSE
##' @param ... optional graphical parameters
plot_ss <- function(time, fit, add = FALSE, ...) {
  ci <- apply(as.matrix(fit$predict), 2, quantile, c(0.025, 
    0.5, 0.975))
  if (!add) {
    args           <- list(...)
    args[["x"]]    <- time
    args[["y"]]    <- ci[2, ]
    args[["type"]] <- "n"
    if (!("ylim" %in% names(args))) 
      args[["ylim"]] <- range(ci, na.rm = TRUE)
    do.call(plot, args)
  }
  ciEnvelope(time, ci[1, ], ci[3, ], col = "lightBlue")
  lines(time, ci[2, ], col = "blue")
  if (!is.null(fit$data)) {
    points(time, fit$data$OBS, pch = "+", cex = 0.5)
  }
}

##' @name solar_geom
##' @title solar_geom
##' @export
##' @author Mike Dietze
##' @description calculates potential top-of-atmosphere shortwave radiation as a function of day of year and location
##' @param doy time as day of year. Integers indicated midnight.
##' @param lon longitude
##' @param lat latitude
solar_geom <- function(doy, lon, lat) {
  
  dt <- median(diff(doy)) * 86400
  hr <- (doy - floor(doy)) * 24
  
  ## calculate potential radiation
  f <- pi/180 * (279.5 + 0.9856 * doy)
  et <- (-104.7 * sin(f) + 596.2 * sin(2 * f) + 4.3 * sin(4 * 
    f) - 429.3 * cos(f) - 2 * cos(2 * f) + 19.3 * cos(3 * 
    f))/3600  #equation of time -> eccentricity and obliquity
  merid <- floor(lon/15) * 15
  merid[merid < 0] <- merid[merid < 0] + 15
  lc <- (lon - merid) * -4/60  ## longitude correction
  tz <- merid/360 * 24  ## time zone
  midbin <- 0.5 * dt/86400 * 24  ## shift calc to middle of bin
  t0 <- 12 + lc - et - tz - midbin  ## solar time
  h <- pi/12 * (hr - t0)  ## solar hour
  dec <- -23.45 * pi/180 * cos(2 * pi * (doy + 10)/365)  ## declination
  
  cosz <- sin(lat * pi/180) * sin(dec) + cos(lat * pi/180) * 
    cos(dec) * cos(h)
  cosz[cosz < 0] <- 0
  
  rpot <- 1366 * cosz
  return(rpot)
}

### NEED TO ADD FUNCTION TO DETECT AND TRIM BURN-IN
