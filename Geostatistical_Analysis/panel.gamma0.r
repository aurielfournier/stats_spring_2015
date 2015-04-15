
panel.gamma0 <- function(x,y,gamma0=gamma0,span=2/3, ...){
  lofit <- loess.smooth(x,y,span=span)
  panel.xyplot(x,y,...)
  panel.xyplot(lofit$x,lofit$y,type="l")
  nump <- round(length(lofit$x)/2)
  dist0 <- approx(lofit$y[1:nump],lofit$x[1:nump],xout=gamma0)$y
  panel.segments(0,gamma0,dist0,gamma0)
  panel.segments(dist0,0,dist0,gamma0)
  parusr <- par()$usr
  panel.text(max(x),min(y),paste("d0=",format(round(dist0,4)),sep=""),adj=1)
}