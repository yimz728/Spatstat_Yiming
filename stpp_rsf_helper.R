## Spatio-temporal point process RSF helper functions


### Sample spatial quadrature points and find tile areas for quadrature weights
###
### Cross spatial and temporal quadratures to obtain 3d locations and volumes
###
SpatTempQuadrature <- function(track.data, environ.data, time.col, time.int, tile.dim){
  require(sp)
  require(rgeos)
  
  environ.data$cellID <- 1:nrow(environ.data)
  # make quad times
  times <- track.data@data[,time.col] - min(track.data@data[,time.col])
  qt <- seq(0, ceiling(max(times)/time.int)*time.int, time.int)
  qt <- sort(c(times, qt))
  qt <- qt[!duplicated(qt)]
  lengths <- diff(c(qt[1], sapply(c(2:length(qt)), function(i,v){mean(v[(i-1):i])}, v=qt), qt[length(qt)]))
  lengths[2] <- sum(lengths[1:2])
  temp.quad <- data.frame(times=qt[-1], length=lengths[-1], obs=ifelse(qt%in%times, 1, 0)[-1])
  temp.quad$x <- temp.quad$y <- NA
  temp.quad$x[temp.quad$obs==1] <- coordinates(track.data)[-1,1]
  temp.quad$y[temp.quad$obs==1] <- coordinates(track.data)[-1,2]
  
  # Get disperion boundries
  delta <- diff(times)
  d1 <- spDists(track.data)
  pwd <- d1[col(d1)==(row(d1)+1)]
  disp.per.time <- pwd/delta
  cut.vls <- quantile(delta, probs=seq(0,1,0.2))
  cut.vls[6] <- cut.vls[6]*1.1
  cut.vls[1] <- 0
  delta.cat <- cut(delta, cut.vls)
  disp.rad.lookup <- aggregate(disp.per.time, list(delta.cat), FUN=max)
  disp.rad.lookup$x <- disp.rad.lookup$x/tile.dim[1]
  disp.rad.lookup$y <- (disp.rad.lookup$x*tile.dim[1])/tile.dim[2]
  
  # Make spatial quad points at temporal quad points
  cat("\nConstructing spatial x temporal quadrature data\n")
  cat("Please be patient...\n")
  mu <- track.data@coords[1,]
  time.last <- 0
  df <- NULL
  quad.time <- temp.quad$time
  time.last <- 0
  pb <- txtProgressBar(min = 0, max = nrow(temp.quad), style = 3)
  for(i in 1:nrow(temp.quad)){
    delta.last <- quad.time[i]-time.last
    rad.x <- ceiling(disp.rad.lookup$x[findInterval(delta.last, cut.vls)]*delta.last)
    rad.y <- ceiling(disp.rad.lookup$y[findInterval(delta.last, cut.vls)]*delta.last)
    x.grid <- c(seq(mu[1]-rad.x*tile.dim[1], mu[1]-tile.dim[1], tile.dim[1]), mu[1], seq(mu[1]+tile.dim[1], mu[1]+rad.x*tile.dim[1], tile.dim[1]))
    y.grid <-  c(seq(mu[2]-rad.y*tile.dim[2], mu[2]-tile.dim[2], tile.dim[2]), mu[2], seq(mu[2]+tile.dim[2], mu[2]+rad.y*tile.dim[2], tile.dim[2]))
    #grd.bbox <- matrix(c(mu[1]-(rad.x+0.5)*tile.dim[1], mu[2]-(rad.y+0.5)*tile.dim[2], mu[1]+(rad.x+0.5)*tile.dim[1], mu[2]+(rad.y+0.5)*tile.dim[2]), 2, 2)
    #colnames(grd.bbox) <- c("min", "max")
    #rownames(grd.bbox) <- c("x","y")
    spqt <- as(SpatialPoints(expand.grid(x=x.grid, y=y.grid), proj4string=CRS(proj4string(environ.data))), "SpatialPixels")
    spqtPoly <- as(spqt, "SpatialPolygons")
    if(temp.quad$obs[i]==1) spqt <- rbind(SpatialPoints(temp.quad[i,c("x","y")], proj4string=CRS(proj4string(environ.data))), as(spqt,"SpatialPoints"))
    else spqt <- as(spqt, "SpatialPoints")
    
    #spqt@bbox <- grd.bbox    
    out <- as(spqt, "data.frame")
    out$t <- temp.quad$time[i]
    out$delta.last <- delta.last
    out$area <- prod(tile.dim)
    out$length<- temp.quad$length[i]
    out$volume <- out$area*temp.quad$length[i]
    out$response <- rep(0,length(out$area))
    if(temp.quad$obs[i]==1){
      out$area <- out$area/c(2, table(spqt%over%spqtPoly))
      out$response[1] <- 1
    } 
    out$response <- out$response/out$volume[1]  
    cov <- spqt %over% environ.data
    out  <- cbind(out, cbind(dist=spDistsN1(spqt, mu), cov))
    out$bmKern <- -0.5*out$dist^2/out$delta.last
    out <- out[!is.na(cov[,1]),]
    df <- rbind(df, out)
    if(temp.quad$obs[i]==1){
      mu <- as.numeric(temp.quad[i,c("x","y")])
      time.last <- temp.quad$time[i]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  class(df) <- c("stQuad", class(df))
  return(df)
}

  


###
### Fit RSF model with STPP
###

stppRSF <- function(model, Q.data, use.gam=FALSE, ...){
  if(!inherits(Q.data, "stQuad")) stop("\nQ.data is not an object of class 'stQuad'! See function 'SpatTempQuadrature'\n")
  if(!use.gam){
    fit <- glm(model, data=Q.data, family = "poisson", weights=volume, ...)
    fit$aic <- 2*(deviance(fit)/2 + sum(log(fit$prior.weights[fit$y!=0])) + sum(fit$y!=0)) + 2*length(fit$coef)
  }
  else{
    if(nrow(Q.data<10000)){
      fit <- gam(model, data=Q.data, family="poisson", weights=volume, ...)
    } else{
      fit <- bam(model, data=Q.data, family="poisson", weights=volume, ...)
    }
    fit$aic <- 2*(deviance(fit)/2 + sum(log(fit$prior.weights[fit$y!=0])) + sum(fit$y!=0)) + 2*sum(fit$edf)
  }
  fit$maxLogLik <- -(deviance(fit)/2 + sum(log(fit$prior.weights[fit$y!=0])) + sum(fit$y!=0))
  class(fit) <- c("stppRSF.fit", class(fit))
  return(fit)
}


summaryRSF <- function(object){
  if(!(inherits(object, "spatRSF.fit") | inherits(object, "stppRSF.fit"))) stop("\nThis is a summery function for 'stppRSF.fit' and 'spatRSF.fit' objects only/n")
  return(summary(object, disp=1))
}













