# The rose function plots a rose diagram with the lengths of the rays given at
# the end of each ray.  This function takes two arguments: a variogram object
# and the value of the semivariogram (gamma) at which the rose plot is desired.
# =============================================================================
rose <- function(dat.var,critgam){
  numcases <- length(dat.var$dist)
  numdirec <- length(unique(dat.var$dir.hor))
  rose.len <- rep(0,numdirec)
  rose.azi <- rep(0,numdirec)
  j <- 1
  for (i in 1:(numcases-1)){
    prod <- (dat.var$gamma[i] - critgam)*(dat.var$gamma[i+1] - critgam)
    if (dat.var$dir.hor[i] != dat.var$dir.hor[i+1]) prod <- 0
    if (j > 1){
      if (rose.azi[j-1] == dat.var$dir.hor[i]) prod <- 0
    }
    if (prod < 0){
        rose.len[j] <- dat.var$dist[i] + (critgam - dat.var$gamma[i])*
                       (dat.var$dist[i+1] - dat.var$dist[i])/
                       (dat.var$gamma[i+1] - dat.var$gamma[i])
        rose.azi[j] <- as.numeric(unique(dat.var$dir.hor)[j])
        j <- j + 1
    }
  }
  rose.azi <- rose.azi*pi/180
  dscale <- max(rose.len)
  rose.len <- rose.len/dscale
  plot(0,0,type="n",axes=F,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),
    xlab="E-W Direction",ylab="",cex.lab=1.6)
  for (j in 1:(length(rose.len))){
    segments(0,0,(rose.len[j]*sin(rose.azi[j])),(rose.len[j]*cos(rose.azi[j])))
    text((rose.len[j]*sin(rose.azi[j]) + (.1*sin(rose.azi[j]))),
         (rose.len[j]*cos(rose.azi[j]) + (.1*cos(rose.azi[j]))),
         (format(round(rose.len[j]*dscale,2))),cex=1.5)
    segments(0,0,-(rose.len[j]*sin(rose.azi[j])),
                 -(rose.len[j]*cos(rose.azi[j])))
  }
  maxx <- max(rose.len*sin(rose.azi))
  text(-(maxx+.25),0,"N-S Direction",srt=90,crt=90,cex=1.6)
  title("Rose Diagram",srt=0,crt=0,cex.main=2.0)
  print("Rose Lengths")
  rose.len*dscale
}
