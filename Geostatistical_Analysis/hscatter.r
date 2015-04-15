hscatter <- function(x,y,u,v,h,dtol=NA,atol=NA){

# This function generates h-scatterplots or cross h-scatterplots in any direction
# desired, for use with lattice and geostatistical data.  Arguments to this
# function are:
#     x = a vector of the x-coordinates of the grid
#     y = a vector of the y-coordinates of the grid
#     u = a vector of the response variable U values
#     v = a vector of the response variable V values
#     h = a 2x1 vector indicating the direction h
#  dtol = a distance tolerance (used with non-lattice data)
#  atol = an angle tolerance (used with non-lattice data)
# The function creates the h-scatterplot and returns the covariance, correlation,
# and semivariogram of the paired values.  If an h-scatterplot is desired, u & v
# should be the same variable, and if a cross h-scatterplot is desired, u & v
# should be different.

 if (missing(dtol) && missing(atol)){
  umat <- tapply(u,list(y,x),function(z)z)
  vmat <- tapply(v,list(y,x),function(z)z)
  if (h[1] > 0) vmat <- vmat[,-((max(x)+1-h[1]):max(x))]
  if (h[1] > 0) umat <- umat[,-(1:h[1])]
  if (h[2] > 0) vmat <- vmat[-((max(y)+1-h[2]):max(y)),]
  if (h[2] > 0) umat <- umat[-(1:h[2]),]
  if (h[1] < 0) vmat <- vmat[,-(1:-h[1])]
  if (h[1] < 0) umat <- umat[,-((max(x)+1+h[1]):max(x))]
  if (h[2] < 0) vmat <- vmat[-(1:-h[2]),]
  if (h[2] < 0) umat <- umat[-((max(y)+1+h[2]):max(y)),]
  umat1 <- tapply(u,list(y,x),function(z)z)
  vmat1 <- tapply(v,list(y,x),function(z)z)
  if (h[1] > 0) vmat1 <- vmat1[,-(1:h[1])]
  if (h[1] > 0) umat1 <- umat1[,-((max(x)+1-h[1]):max(x))]
  if (h[2] > 0) vmat1 <- vmat1[-(1:h[2]),]
  if (h[2] > 0) umat1 <- umat1[-((max(y)+1-h[2]):max(y)),]
  if (h[1] < 0) vmat1 <- vmat1[,-((max(x)+1+h[1]):max(x))]
  if (h[1] < 0) umat1 <- umat1[,-(1:-h[1])]
  if (h[2] < 0) vmat1 <- vmat1[-((max(y)+1+h[2]):max(y)),]
  if (h[2] < 0) umat1 <- umat1[-(1:-h[2]),]
  semivar <- (1/(2*length(c(vmat))))*sum((c(umat) -
             c(umat1))*(c(vmat1) - c(vmat)))
 }
 else{
  n <- length(u)
  dmat <- matrix(0,nrow=n,ncol=n)
  xy <- cbind(x,y)
  dmat[lower.tri(dmat)] <- dist(xy)
  dmat <- dmat + t(dmat)
  xdiff <- matrix(rep(x,n),nrow=n,byrow=T) - matrix(rep(x,n),nrow=n)
  ydiff <- matrix(rep(y,n),nrow=n,byrow=T) - matrix(rep(y,n),nrow=n)
  angle <- atan(ydiff/xdiff)*180/pi
  disth <- sqrt(sum(h*h))
  angleh <- atan(h[2]/h[1])*180/pi
  index <- seq(1,n*n)
  dmat <- c(dmat); angle <- c(angle)
  iloc <- index[(abs(dmat-disth)<dtol)&(abs(angle-angleh)<atol)]
  irow <- trunc(((iloc-1)/n)+1)
  icol <- iloc - trunc((iloc-1)/n)*n
  if (h[2]==0){
    keep <- irow<icol
    irow <- irow[keep]; icol <- icol[keep]
  }
  vmat <- v[icol]; umat <- u[irow]
  semivar <- (1/(2*length(vmat)))*sum((u[irow]-u[icol])*(v[irow]-v[icol]))
 }
 xlims <- c(min(u),max(u))
 ylims <- c(min(v),max(v))
 lab1 <- paste(deparse(substitute(u)),"(t)",sep="")
 lab2 <- paste(deparse(substitute(v)),"(t+h)",sep="")
 plot(umat,vmat,pch="+",xlab="",ylab="",xlim=xlims,ylim=ylims,axes=F,cex=1.2)
 axis(1,tck=.02,pos=min(u,v),mgp=c(3,.5,0),cex.axis=1.3)
 axis(2,tck=.02,pos=min(u,v),mgp=c(3,.5,0),cex.axis=1.3)
 mtext(side=1,line=2.0,lab1,cex=1.3)
 mtext(side=2,line=2.3,lab2,cex=1.3)
 mtext(side=3,line=0,paste("h = (",h[1],",",h[2],")",sep=""),cex=1.3)
 lines(c(min(u,v),max(u,v)),c(min(u,v),max(u,v)),lty=8)
 numpairs <- length(c(vmat))
 cat("\n Correlation:   p(h)=", cor(c(vmat),c(umat)))
 cat("\n Covariance:    C(h)=", cor(c(vmat),c(umat))*sqrt(var(c(umat))*
                              var(c(vmat)))*(numpairs-1)/numpairs)
 cat("\n Semivariogram: g(h)=", semivar)
 cat("\n")
 data.frame(covar=cor(c(vmat),c(umat))*sqrt(var(c(umat))*var(c(vmat)))*
   (numpairs-1)/numpairs,corr=cor(c(vmat),c(umat)),semivar=semivar)
}

