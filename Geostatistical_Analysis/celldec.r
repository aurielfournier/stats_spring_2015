# The function cell.dec performs cell declustering to estimate the
# global mean over some region.  There are 5 arguments to this
# function, as outlined below, and the function returns the global
# mean estimate.  The five arguments are:
#       x = the vector of x-coordinates of sample values,
#       y = the vector of y-coordinates of sample values,
#       v = the vector of response values at the sample locations,
#  xcells = the number of cells desired in the x-direction,
#  ycells = the number of cells desired in the y-direction.
#
# As an example, if we wanted to perform cell declustering on the
# Walker Lake data with 3 cells in the x-direction, and 5 cells in
# the y-direction, we would type:
#    > attach(walk470)
#    > globmean <- celldec(x,y,v,3,5)

celldec <- function(x,y,v,xcells,ycells){
  leftlim <- min(x)-(x[order(x)[2]]-min(x))/2
  rightlim <- max(x)+(max(x)-x[order(x)[length(x)-1]])/2
  downlim <- min(y)-(y[order(y)[2]]-min(y))/2
  uplim <- max(y)+(max(y)-y[order(y)[length(y)-1]])/2
  xgap <- (rightlim-leftlim)/xcells
  ygap <- (uplim-downlim)/ycells
  xmids <- (leftlim+rightlim)/2; ymids <- (downlim+uplim)/2
  if (xcells!=1) xmids <- seq(leftlim+xgap/2,rightlim-xgap/2,xgap)
  if (ycells!=1) ymids <- seq(downlim+ygap/2,uplim-ygap/2,ygap)
  xmat <- t(matrix(rep(x,length(xmids)),nrow=length(x)))
  ymat <- t(matrix(rep(y,length(ymids)),nrow=length(y)))
  xmids <- matrix(rep(xmids,length(x)),nrow=length(xmids))
  ymids <- matrix(rep(ymids,length(y)),nrow=length(ymids))
  distx <- matrix(abs(xmat-xmids),nrow=xcells)
  disty <- matrix(abs(ymat-ymids),nrow=ycells)
  xcelllab <- matrix(apply(distx,2,order),nrow=xcells)[1,]
  ycelllab <- matrix(apply(disty,2,order),nrow=ycells)[1,]
  celllab <- (ycelllab-1)*xcells + xcelllab
  d <- length(unique(celllab))
  numrep <- rep(0,length(x))
  for (i in 1:length(x)) numrep[i] <- sum(celllab==celllab[i])
  cdmean <- sum((1/d)*(v/numrep))
}
