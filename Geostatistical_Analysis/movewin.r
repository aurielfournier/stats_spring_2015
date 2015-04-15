movewin <- function(x=x,y=y,v=v,wx=wx,wy=wy,wo=wo){

# This function computes moving window means and variances for any type of data.
# Arguments to this function are:
#      x = a vector of the x-coordinates of the data
#      y = a vector of the y-coordinates of the data
#      v = a vector of the response variable values
#     wx = width of the window in the x-direction
#     wy = width of the window in the y-direction
#     wo = width of the overlap between adjacent windows
# The function returns a data frame with new (x,y) coordinates corresponding to the
# midpoints of the windows, the number of points in each window (numvals), the
# moving window means (means), and the moving window standard deviations (sdevs).

  edge <- min(x[x!=min(x)]-min(x))/2
  dimx <- max(x) - min(x) + 2*edge
  dimy <- max(y) - min(y) + 2*edge
  numx <- trunc((dimx-wx)/(wx-wo) + 1)
  numy <- trunc((dimy-wy)/(wy-wo) + 1)
  xmid <- matrix(rep(min(x)-edge+(0:(numx-1))*(wx-wo),numy),nrow=numx)+(wx/2)
  ymid <- matrix(rep(min(y)-edge+(0:(numy-1))*(wy-wo),numx),nrow=numx,byrow=T)+(wy/2)
  isin <- array(dim=c(numx,numy,length(x)))
  for (i in 1:length(x)) isin[,,i] <- (abs(x[i]-xmid)<=(wx/2))+(abs(y[i]-ymid)<=(wy/2))==2
  numvals <- apply(isin,c(1,2),sum)
  means <- matrix(nrow=numx,ncol=numy)
  sdevs <- matrix(nrow=numx,ncol=numy)
  for (i in 1:numx){
    for (j in 1:numy){
      if (numvals[i,j]>0){
        means[i,j] <- mean(v[isin[i,j,]==T])
        sdevs[i,j] <- sd(v[isin[i,j,]==T])
      }
    }
  }
  data.frame(x=c(xmid),y=c(ymid),numvals=c(numvals),means=c(means),sdevs=c(sdevs))
}
