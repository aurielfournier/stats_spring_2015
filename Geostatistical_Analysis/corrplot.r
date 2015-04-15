corrplot <- function(x,y,u,v,h,numlags,dtol=NA,atol=NA){

# This function generates the covariogram, correlogram, and semivariogram plots in the
# direction specified by "h", and with "numlags" lags in that direction.  Additionally,
# the values of these functions are returned as a data frame.  Arguments to this
# function are:
#         x = a vector of the x-coordinates of the grid
#         y = a vector of the y-coordinates of the grid
#         u = a vector of the response variable U values
#         v = a vector of the response variable V values
#         h = a 2x1 vector indicating the direction h
#   numlags = the number of lags at which to determine the functions
#      dtol = a distance tolerance for non-lattice data (optional)
#      atol = an angle tolerance for non-lattice data (optional)
# If the cross-covariogram, cross-correlogram, and cross-semivariogram are desired,
# then u & v should be different.  Otherwise, u & v should be the same.

correlog <- rep(0,numlags)
covariog <- rep(0,numlags)
semivar <- rep(0,numlags)

for (i in 1:numlags){
  h2 <- h*i
  if (missing(dtol)) dtol <- max(abs(h2))*.0001
  if (missing(atol)) atol <- .0001
  n <- length(u)
  dmat <- matrix(0,nrow=n,ncol=n)
  xy <- cbind(x,y)
  dmat[lower.tri(dmat)] <- dist(xy)
  dmat <- dmat + t(dmat)
  xdiff <- matrix(rep(x,n),nrow=n,byrow=T) - matrix(rep(x,n),nrow=n)
  ydiff <- matrix(rep(y,n),nrow=n,byrow=T) - matrix(rep(y,n),nrow=n)
  angle <- atan(ydiff/xdiff)*180/pi
  disth <- sqrt(sum(h2*h2))
  angleh <- atan(h2[2]/h2[1])*180/pi
  index <- seq(1,n*n)
  dmat <- c(dmat); angle <- c(angle)
  iloc <- index[(abs(dmat-disth)<=dtol)&(abs(angle-angleh)<=atol)]
  irow <- trunc(((iloc-1)/n)+1)
  icol <- iloc - trunc((iloc-1)/n)*n
  if (h[2]==0){
    keep <- irow<icol
    irow <- irow[keep]; icol <- icol[keep]
  }
  vmat <- v[icol]; umat <- u[irow]
  nh <- sum(!is.na((u[irow]-u[icol])*(v[irow]-v[icol])))
  semivar[i] <- (1/(2*nh))*sum((u[irow]-u[icol])*(v[irow]-v[icol]),na.rm=T)
  numpairs <- min(sum(!is.na(vmat)),sum(!is.na(umat)))
  umatvec <- c(umat); vmatvec <- c(vmat)
  include <- !is.na(vmatvec)&!is.na(umatvec)
  correlog[i] <- cor(vmatvec[include],umatvec[include])
  covariog[i] <- cor(vmatvec[include],umatvec[include])*
    sqrt(var(umatvec,na.rm=T)*var(vmatvec,na.rm=T))*(numpairs-1)/numpairs
}
include <- !is.na(v)&!is.na(u)
covar0 <- cor(v[include],u[include])*sqrt(var(u,na.rm=T)*var(v,na.rm=T))*
  (numpairs-1)/numpairs
corr0 <- cor(v[include],u[include])
if (h[1]==0 && h[2]>0) azimuth <- 0
if (h[1]==0 && h[2]<0) azimuth <- 180
if (h[1]!=0) azimuth <- 90-atan(h[2]/h[1])*180/pi

par(mfrow=c(2,2))
xlims <- c(0,numlags); ylims <- c(min(covariog,covar0,0),max(covariog,covar0))
plot(0:numlags,c(covar0,covariog),type="b",xlab="",ylab="",xlim=xlims,ylim=ylims,
  axes=F,cex=1.0)
axis(1,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
axis(2,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
mtext(side=1,line=2.0,"Lag (h)",cex=1.0)
mtext(side=2,line=2.0,"Covariogram C(h)",cex=1.0)
xlims <- c(0,numlags); ylims <- c(min(correlog,corr0,0),max(correlog,corr0))
plot(0:numlags,c(corr0,correlog),type="b",xlab="",ylab="",xlim=xlims,ylim=ylims,
  axes=F,cex=1.0)
axis(1,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
axis(2,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
mtext(side=1,line=2.0,"Lag (h)",cex=1.0)
mtext(side=2,line=2.7,"Correlogram p(h)",cex=1.0)
xlims <- c(0,numlags); ylims <- c(0,max(semivar))
plot(0:numlags,c(0,semivar),type="b",xlab="",ylab="",xlim=xlims,ylim=ylims,
  axes=F,cex=1.0)
axis(1,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
axis(2,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
mtext(side=1,line=2.0,"Lag (h)",cex=1.0)
mtext(side=2,line=2.7,"Semivariogram Gamma(h)",cex=1.0)
plot(c(0,1),c(0,1),xlab="",ylab="",axes=F,pch=" ")
text(.5,.6,paste("Azimuth = ",round(azimuth,1),sep=""),cex=1.5)

data.frame(lag=1:numlags,corr=correlog,covar=covariog,semivar=semivar)
}

