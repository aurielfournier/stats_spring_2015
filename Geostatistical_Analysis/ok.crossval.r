ok.crossval <- function(x, y, v, nugget, sill, range, azim, ratio, 
                        model,r,p){

# This function performs crossvalidation for ordinary kriging, polygonal
# declustering, triangulation, inverse distance weighting, and the local
# sample mean estimation methods.  Arguments to this function are:
#     x      = a vector of the x-coordinates of the locations
#     y      = a vector of the y-coordinates of the locations
#     v      = a vector of the response values at the locations
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
#         r  = the radius of influence for local sample mean
#         p  = the power of the inverse weighting
#
# If the variogram model is isotropic, choose the range ratio to be 0.
#
# This function requires the libraries "tripack" and "splancs" and uses
#   the "ordkrige" and "getCD" functions.
# The function prints out several crossvalidation statistics for each estimation
# method and returns five vectors of the crossvalidation values for the five
# estimation procedures.  They are stored in a data frame with field names
# pred.pd, pred.tr, pred.ls, pred.iw, and pred.ok.

  n <- length(x)
  pred.ok <- rep(0,n); pred.ls <- rep(0,n)
  pred.iw <- rep(0,n); pred.tr <- rep(NA,n)
  dist.all <- dist(cbind(x,y))
  dist.mat <- matrix(0,attr(dist.all,"Size"),attr(dist.all,"Size"))
  dist.mat[lower.tri(dist.mat)] <- dist.all
  dist.mat <- dist.mat + t(dist.mat)
  pred.pd <- apply(dist.mat,1,function(a,v) v[order(a)[2]],v)
  covs <- getCD(x,y,x,y,nugget,sill,range,azim,ratio,model)
  Cinv <- solve(covs$C)
  for(i in 1:n){
    cat("i = ",i,"\n")
    dist0 <- dist.mat[i,]
    local <- seq(1,length(dist0))[-c(seq(1,length(dist0))[dist0>r],i)]
    if (length(local)<1) stop("Radius of Influence Too Small")
    pred.ls[i] <- mean(v[local])
    weights <- (1/dist0[local]**p)/sum(1/dist0[local]**p)
    pred.iw[i] <- sum(v[local]*weights)
    xi <- x[-i]; yi <- y[-i]; vi <- v[-i]
    if (!is.element(i,chull(x,y))){
      delaun <- triangles(tri.mesh(xi,yi))
      numtriang <- length(delaun[,1])
      tri <- as.vector(tri.find(tri.mesh(xi,yi),x[i],y[i]),mode="integer")
      thistri <- rep(0,numtriang)
      for (j in 1:numtriang) thistri[j] <- setequal(tri,delaun[j,1:3])
      trinum <- rev(order(thistri))[1]
      A1 <- areapl(cbind(c(xi[tri[1:2]],x[i]),c(yi[tri[1:2]],y[i])))
      A2 <- areapl(cbind(c(xi[tri[c(1,3)]],x[i]),c(yi[tri[c(1,3)]],y[i])))
      A3 <- areapl(cbind(c(xi[tri[2:3]],x[i]),c(yi[tri[2:3]],y[i])))
      A <- A1+A2+A3
      mina <- min(A1,A2,A3)/A
      if ((!is.element(0,delaun[trinum,4:6]))&&(mina>.025))
        pred.tr[i] <- (A1*vi[tri[3]] + A2*vi[tri[2]] + A3*vi[tri[1]])/A
    }
    krige.out <- ordkrige1(v[ - i], Cinv, covs$D[,i], sill,i)
    pred.ok[i] <- krige.out$pred
  }
  pred.ok[pred.ok<min(v)] <- min(v)
  pred.ok[pred.ok>max(v)] <- max(v)
  nt <- n - sum(is.na(pred.tr))
  nonna <- order(pred.tr)[1:nt]
  mse.pd <- (1/n)*sum((pred.pd - v)^2)
  mse.ls <- (1/n)*sum((pred.ls - v)^2)
  mse.iw <- (1/n)*sum((pred.iw - v)^2)
  mse.tr <- (1/nt)*sum((pred.tr[nonna]-v[nonna])**2)
  mse.ok <- (1/n)*sum((pred.ok - v)^2)
  mae.pd <- (1/n)*sum(abs(pred.pd - v))
  mae.ls <- (1/n)*sum(abs(pred.ls-v))
  mae.iw <- (1/n)*sum(abs(pred.iw-v))
  mae.tr <- (1/nt)*sum(abs(pred.tr[nonna]-v[nonna]))
  mae.ok <- (1/n)*sum(abs(pred.ok - v))
  corr.pd <- cor(pred.pd,v)
  corr.ls <- cor(pred.ls,v)
  corr.iw <- cor(pred.iw,v)
  if (nt>1) corr.tr <- cor(pred.tr[nonna],v[nonna])
  if (nt<=1) corr.tr <- NA
  corr.ok <- cor(pred.ok,v)
  numdec <- trunc(log10(mae.pd) - 1)
  cat("Polygonal Declustering","\n======================\n")
  cat("n = ", n, "   True Mean = ", round(mean(v),numdec),"   m =",
      round(mean(pred.pd),numdec), "   SD = ", round(sqrt(var(
      pred.pd)), numdec), "\n")
  cat("5-Number Summary (True Values): (", round(quantile(v),
      numdec), ")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(
      pred.pd), numdec), ")\n")
  cat("5-Number Summary (Residuals): (", round(
      quantile(pred.pd - v), numdec), ")\n")
  cat("Residual (Mean,SD) = (", round(mean(pred.pd - v), numdec),
      round(sqrt(var(pred.pd - v)), numdec), ")\n")
  cat("MSE = ", round(mse.pd, numdec),"   MAE = ", round(mae.pd,
      numdec),"   Correlation = ", round(corr.pd, 3), "\n\n")
  cat("Delaunay Triangulation","\n======================\n")
  cat("n = ",nt,"   True Mean = ",round(mean(v[nonna]),numdec),"   m = ",round(mean
    (pred.tr[nonna]),numdec),"   SD = ",round(sqrt(var(pred.tr[nonna])),
    numdec),"\n")
  cat("5-Number Summary (True Values): (",round(quantile(v[nonna]),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.tr[nonna]),
    numdec),")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.tr[nonna]-
    v[nonna]),numdec),")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.tr[nonna]-v[nonna]),numdec),
    round(sqrt(var(pred.tr[nonna]-v[nonna])),numdec),")\n")
  cat("MSE = ",round(mse.tr,numdec),"   MAE = ",round(mae.tr,numdec),
      "   Correlation = ",round(corr.tr,3),"\n\n")
  cat("Local Sample Mean (Unweighted)","\n==============================\n")
  cat("n = ",n,"   True Mean = ",round(mean(v),numdec),"   m = ",
    round(mean(pred.ls),numdec),"   SD = ",round(sqrt(var(pred.ls)),numdec),
    "   p = ",p,"   r = ",r,"\n")
  cat("5-Number Summary (True Values): (",round(quantile(v),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.ls),numdec),
    ")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.ls-v),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.ls-v),numdec),
    round(sqrt(var(pred.ls-v)),numdec),")\n")
  cat("MSE = ",round(mse.ls,numdec),"   MAE = ",round(mae.ls,numdec),
      "   Correlation = ",round(corr.ls,3),"\n\n")
  cat("Inverse Distance Weighting","\n==========================\n")
  cat("n = ",n,"   True Mean = ",round(mean(v),numdec),"   m = ",
    round(mean(pred.iw),numdec),"   SD = ",round(sqrt(var(pred.iw)),numdec),
    "   p = ",p,"   r = ",r,"\n")
  cat("5-Number Summary (True Values): (",round(quantile(v),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.iw),numdec),
    ")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.iw-v),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.iw-v),numdec),
    round(sqrt(var(pred.iw-v)),numdec),")\n")
  cat("MSE = ",round(mse.iw,(numdec+1)),"   MAE = ",round(mae.iw,(numdec+1)),
      "   Correlation = ",round(corr.iw,3),"\n\n")
  cat("Ordinary Kriging", "\n======================\n")
  cat("n = ", n, "   True Mean = ", round(mean(v),numdec), "   m = ",
      round(mean(pred.ok),numdec), "   SD = ", round(sqrt(var(
      pred.ok)), numdec), "\n")
  cat("5-Number Summary (True Values): (", round(quantile(v),
      numdec), ")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(
      pred.ok), numdec), ")\n")
  cat("5-Number Summary (Residuals): (", round(quantile(pred.ok - v),
      numdec), ")\n")
  cat("Residual (Mean,SD) = (", round(mean(pred.ok - v), numdec),
      round(sqrt(var(pred.ok - v)), numdec), ")\n")
  cat("MSE = ", round(mse.ok, numdec),"   MAE = ", round(mae.ok,
      numdec),"   Correlation = ", round(corr.ok, 3), "\n\n")
  list(pred.pd=pred.pd,pred.tr=pred.tr,pred.ls=pred.ls,pred.iw=pred.iw,
       pred.ok=pred.ok)
}
