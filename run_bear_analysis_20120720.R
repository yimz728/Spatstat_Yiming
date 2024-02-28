# Bear data analysis
# written by: Devin Johnson, email: devin.johnson@noaa.gov
library(sp)
library(mgcv)
source("stpp_rsf_helper.R")

habmat = read.table("bear_habitat.csv",head=T,sep=',')
habmat$x <- habmat$x-0.5
habmat$y <- habmat$y-0.5
coordinates(habmat) <- ~ x + y
habmat <- as(habmat, "SpatialPixelsDataFrame")

track <- read.table("bear_track.csv",head=T,sep=',')
coordinates(track) <- ~x+y

 # Obtain spatial quadrature points and tile areas
stQuadData <- SpatTempQuadrature(track.data=track, environ.data=habmat, time.col="time", time.int=10, tile.dim=c(1,1))

fit <- stppRSF(response ~  strDist + bmKern + s(x,y, k=25, fx=TRUE), Q.data=stQuadData, use.gam=TRUE)
summ <- summaryRSF(fit)



save(list=ls(), file="bear_stpp_results.RData")
  




