# The function point.crossval takes as arguments a vector of x-coordinates (x),
# a vector of y-coordinates (y), a vector of response values (z), a power for
# inverse distance weighting (p) and a radius of influence for both inverse
# distance weighting and the unweighted local sample mean (r).  It performs
# crossvalidation on the observed values using the point estimation methods:
# polygonal declustering, triangulation, unweighted local sample mean, and
# inverse distance weighting.  Tables of output for each of the four methods
# are produced with means, standard deviations, 5-number summaries, MSE, MAE,
# and the correlation between the true and estimated values.  Also returned
# with the function under the names pred.pd, pred.tr, pred.ls, pred.iw are the
# vectors of predicted values corresponding to the true values on which
# crossvalidation was performed.

# Libraries requried:  tripack, splancs

point.crossval <- function(x,y,z,p,r){
  n <- length(x); pred.pd <- rep(0,n); pred.iw <- rep(0,n)
  pred.tr <- rep(NA,n); pred.ls <- rep(0,n)
  for (i in 1:n){
    cat("i = ",i,"\n")
    dist0 <- sqrt((x-x[i])**2 + (y-y[i])**2)
    pred.pd[i] <- z[order(dist0)[2]]
    local <- seq(1,length(dist0))[-c(seq(1,length(dist0))[dist0>r],i)]
    if (length(local)<1) stop("Radius of Influence Too Small")
    pred.ls[i] <- mean(z[local])
    weights <- (1/dist0[local]**p)/sum(1/dist0[local]**p)
    pred.iw[i] <- sum(z[local]*weights)
    if (!is.element(i,chull(x,y))){
      delaun <- triangles(tri.mesh(x[-i],y[-i]))
      numtriang <- length(delaun[,1])
      area <- rep(0,numtriang)
      x2 <- x[-i]; y2 <- y[-i]; z2 <- z[-i]
      intri <- rep(F,numtriang)
      for (j in 1:numtriang)
        intri[j] <- in.convex.hull(tri.mesh(x2[delaun[j,1:3]],
                                   y2[delaun[j,1:3]]),x[i],y[i])
      tri <- rev(order(intri==T))[1]
      A1 <- areapl(cbind(c(x2[delaun[tri,1:2]],x[i]),c(y2[delaun[tri,1:2]],y[i])))
      A2 <- areapl(cbind(c(x2[delaun[tri,c(1,3)]],x[i]),c(y2[delaun[tri,c(1,3)]],y[i])))
      A3 <- areapl(cbind(c(x2[delaun[tri,2:3]],x[i]),c(y2[delaun[tri,2:3]],y[i])))
      A <- A1+A2+A3
      pred.tr[i] <- (A1*z[delaun[tri,3]] + A2*z[delaun[tri,2]] + A3*z[delaun[tri,1]])/A
    }
  }
  nt <- n - sum(is.na(pred.tr))
  nonna <- order(pred.tr)[1:nt]
  mse.pd <- (1/n)*sum((pred.pd-z)**2)
  mse.ls <- (1/n)*sum((pred.ls-z)**2)
  mse.iw <- (1/n)*sum((pred.iw-z)**2)
  mse.tr <- (1/nt)*sum((pred.tr[nonna]-z[nonna])**2)
  mae.pd <- (1/n)*sum(abs(pred.pd-z))
  mae.ls <- (1/n)*sum(abs(pred.ls-z))
  mae.iw <- (1/n)*sum(abs(pred.iw-z))
  mae.tr <- (1/nt)*sum(abs(pred.tr[nonna]-z[nonna]))
  corr.pd <- cor(pred.pd,z)
  corr.ls <- cor(pred.ls,z)
  corr.iw <- cor(pred.iw,z)
  if (nt>1) corr.tr <- cor(pred.tr[nonna],z[nonna])
  if (nt<=1) corr.tr <- NA
  numdec <- trunc(log10(mae.pd)-1)
  cat("Polygonal Declustering","\n======================\n")
  cat("n = ",n,"   True Mean = ",round(mean(z),numdec),"   m = ",round(mean(
    pred.pd),numdec),"   SD = ",round(sqrt(var(pred.pd)),numdec),"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.pd),
    numdec),")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.pd-z),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.pd-z),numdec),
    round(sqrt(var(pred.pd-z)),numdec),")\n")
  cat("MSE = ",round(mse.pd,numdec),"   MAE = ",round(mae.pd,numdec),
      "   Correlation = ",round(corr.pd,3),"\n\n")
  cat("Delaunay Triangulation","\n======================\n")
  cat("n = ",nt,"   True Mean = ",round(mean(z[nonna]),numdec),"   m = ",round(mean
    (pred.tr[nonna]),numdec),"   SD = ",round(sqrt(var(pred.tr[nonna])),
    numdec),"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.tr[nonna]),
    numdec),")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.tr[nonna]-
    z[nonna]),numdec),")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.tr[nonna]-z[nonna]),numdec),
    round(sqrt(var(pred.tr[nonna]-z[nonna])),numdec),")\n")
  cat("MSE = ",round(mse.tr,numdec),"   MAE = ",round(mae.tr,numdec),
      "   Correlation = ",round(corr.tr,3),"\n\n")
  cat("Local Sample Mean (Unweighted)","\n==============================\n")
  cat("n = ",n,"   True Mean = ",round(mean(z),numdec),"   m = ",
    round(mean(pred.ls),numdec),"   SD = ",round(sqrt(var(pred.ls)),numdec),
    "   p = ",p,"   r = ",r,"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.ls),numdec),
    ")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.ls-z),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.ls-z),numdec),
    round(sqrt(var(pred.ls-z)),numdec),")\n")
  cat("MSE = ",round(mse.ls,numdec),"   MAE = ",round(mae.ls,numdec),
      "   Correlation = ",round(corr.ls,3),"\n\n")
  cat("Inverse Distance Weighting","\n==========================\n")
  cat("n = ",n,"   True Mean = ",round(mean(z),numdec),"   m = ",
    round(mean(pred.iw),numdec),"   SD = ",round(sqrt(var(pred.iw)),numdec),
    "   p = ",p,"   r = ",r,"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.iw),numdec),
    ")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.iw-z),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.iw-z),numdec),
    round(sqrt(var(pred.iw-z)),numdec),")\n")
  cat("MSE = ",round(mse.iw,(numdec+1)),"   MAE = ",round(mae.iw,(numdec+1)),
      "   Correlation = ",round(corr.iw,3),"\n\n")
  list(pred.pd=pred.pd,pred.tr=pred.tr,pred.ls=pred.ls,pred.iw=pred.iw)
}
