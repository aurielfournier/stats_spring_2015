getCD <- function(x,y,x0,y0,nugget,sill,range,azim,ratio,model){

# This function finds the matrix of covariances C between all sample locations
# given in the location vectors x,y, and the matrix of covariances D between
# the vector of sample locations (x,y) and the vector of prediction sites
# (x0,y0).  Arguments to this function are:
#     x      = a vector of the x-coordinates of the locations
#     y      = a vector of the y-coordinates of the locations
#     x0     = a vector of the x-coordinates of the k prediction points
#     y0     = a vector of the y-coordinates of the k prediction points
#     nugget = the nugget for the variogram model
#     sill   = the sill for the variogram model (sigma^2)
#     range  = the range for the variogram model
#     azim   = the direction (from 90) of maximum spatial continuity for
#              geometric anisotropy.
#     ratio  = the range ratio associated with the azimuth for geometric
#              anisotropy.
#     model  = the type of variogram model.  Choices are:
#                 "spher" = Spherical
#                 "expon" = Exponential
#                 "gauss" = Gaussian
#
# If the variogram model is isotropic, choose the range ratio to be 0.
#
# The function returns a list with the (n+1)x(n+1) matrix of covariances
# among the sample values labeled C, and the (n+1)xk matrix of covariances
# between the sample values and the prediction points labeled D.  Both
# matrices have the supplemented ones columns supplemented for use in the
# kriging equations.  Each column of values in the D matrix corresponds to
# a different prediction location.

  n <- length(x)
  k <- length(x0)
  xy <- matrix(c(x,x0,y,y0),nrow=(n+k),byrow=F)
  if (ratio > 0){
    az <- (90-azim)*pi/180
    rot <- matrix(c(cos(az),sin(az),-sin(az),cos(az)),nrow=2,byrow=T)
    newxy <- rot%*%t(xy)
    disx <- dist(newxy[1,],method="euclidean")
    disy <- dist(newxy[2,],method="euclidean")
    disxy <- sqrt((disx/range)**2 + (disy*ratio/range)**2)
  }
  else{
    disxy <- dist(xy,method="euclidean")/range
  }
  dmat <- matrix(0,(n+k),(n+k))
  dmat[lower.tri(dmat)] <- disxy
  dmat <- dmat + t(dmat)
  if (model == "spher"){
    vmat <- (nugget*(dmat > 0)) + (sill-nugget)*((1.5*dmat) - (0.5*(dmat)**3))*
            (dmat < 1) + (sill-nugget)*(dmat >= 1)
  }
  if (model == "expon"){
    vmat <- (nugget*(dmat > 0)) + (sill-nugget)*(1 - exp(-3*dmat))
  }
  if (model == "gauss"){
    vmat <- (nugget*(dmat > 0)) + (sill-nugget)*(1 - exp(-3*(dmat**2)))
  }
  
  vmat <- sill - vmat
  C <- rbind(cbind(vmat[1:n,1:n],rep(1,n)),t(rep(1,(n+1))))
  C[n+1,n+1] <- 0
  D <- rbind(matrix(vmat[1:n,(n+1):(n+k)],n,k),t(rep(1,k)))
  list(C=C, D=D)
}
