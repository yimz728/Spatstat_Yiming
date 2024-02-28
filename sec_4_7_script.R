##########################
### Bear STPP analysis...
##########################

library(sp)
library(mgcv)
source("gray.colors.rev.R")
source("draw.contour.R")
source("draw.contour.2.R")
load("bear_stpp_results.RData")


pred_dat2 = as.data.frame(habmat)
pred_dat2$pred = habmat$strDist*coef(fit)[2]
strmpred = matrix(pred_dat2$pred, 256, 256)
strmpred = ifelse(strmpred < -10, -10, strmpred)

# Include the following parts only for computing "breaks"
pred_dat = data.frame(strDist=0, delta.last=1, bmKern=0, coordinates(habmat))
pred_dat$smooth = predict(fit, newdata = pred_dat) - coef(fit)[1] 
smooth = matrix(pred_dat$smooth, ncol=256, nrow=256)
smooth = ifelse(smooth < -10, -10, smooth)
xy = coordinates(habmat)

obsid=101 
obs = track[obsid,] 
delta=track$time[obsid+1]-track$time[obsid] 
pred_dat$bmKern = -coef(fit)[3]*((xy[,1]-obs$x)^2 + (xy[,2]-obs$y)^2)/(2*delta)
kern = matrix(pred_dat$bmKern, 256, 256) 
kern = ifelse(kern < -10, -10, kern)
trans = matrix(pred_dat$bmKern + pred_dat$smooth, 256, 256)
trans = ifelse(trans < -10, -10, trans)
# "smooth" and "trans" don't exist in original article, but they're needed for 
# computing the "breaks" which is an argument by plotting fitted selection surface

breaks = quantile(unique(c(smooth, strmpred, trans)), seq(0,1,length=101))

pdf("bear_stpp.pdf", width=8, height=8)
par(mar=c(1,1,2,1)+.1)
layout(matrix(1:4, ncol=2, byrow=T))
image(matrix(habmat@data$strDist, 256, 256)[,1:256], 
      x=0.5:255.5, y=0.5:255.5,cex.axis=1.2,
      cex.lab=1.5, xlab=NULL, ylab=NULL,
      col = gray.colors.rev(100, 0.1), main="distance to the stream and locations",
      asp=1,xaxt="n",yaxt="n",cex.axis=1.2,cex.lab=1.5,cex.main=1.5, bty="n"
      )
points(track$x,track$y, pch=20, col=rgb(0,0,0,0.4), cex=1.5)

image(strmpred[,1:256], 
      x=0.5:255.5, y=0.5:255.5,cex.axis=1.2,
      cex.lab=1.5, xlab=NULL, ylab=NULL,
      col = gray.colors.rev(100, 0.1), main="fitted selection surface",
      asp=1,xaxt="n",yaxt="n",cex.axis=1.2,cex.lab=1.5,cex.main=1.5, bty="n",
      breaks=breaks
      )
dev.off()
