# Defines a function to create normal PP-plots which takes four arguments:
#        x = the vector of data values,
#   points = the values at which normal quantiles are calculated,
#   xticks = the x-axis tick mark locations desired,
#   xlabel = the x-axis label.
# ========================================================================
ppnorm <- function(x,points,xticks,xlabel){
  repx <- matrix(rep(x,length(points)),nrow=length(points),byrow=T)
  reppt <- matrix(rep(points,length(x)),nrow=length(points),byrow=F)
  diff <- reppt-repx
  cumfreq <- apply(diff,1,function(a) mean(a>=0))
  zq <- qnorm(cumfreq,0,1)
  minx <- 2*points[1]-points[2]
  plot(points,zq,xlab="",ylab="",xlim=c(minx,
    2*points[length(points)]-points[length(points)-1]),ylim=c(-2.5,2.5),
    pch="+",axes=F,cex=1.5,type="b")
  axis(1,at=xticks,tck=.02,pos=-2.5,mgp=c(3,.5,0),cex.axis=1.5)
  axis(2,at=qnorm(c(1,2,5,10,20,50,80,90,95,98,99)/100,0,1),tck=.02,
    pos=minx,mgp=c(3,.5,0),labels=c("1","2","5","10","20","50","80","90",
    "95","98","99"),las=2,cex.axis=1.3)
  mtext(side=1,line=2.3,xlabel,cex=1.5)
  mtext(side=2,line=2.0,"Cumulative Frequency (%)",cex=1.5)
  i25 <- sum(cumfreq<.25)
  end25 <- points[i25] + (points[i25+1]-points[i25])*((.25-cumfreq[i25])/
    (cumfreq[i25+1]-cumfreq[i25]))
  lines(c(minx,end25),c(qnorm(.25,0,1),qnorm(.25,0,1)),lty=8)
  lines(c(end25,end25),c(-2.5,qnorm(.25,0,1)),lty=8)
  text(end25,-2.62,"Q1",cex=1.2)
  i75 <- sum(cumfreq<.75)
  end75 <- points[i75] + (points[i75+1]-points[i75])*((.75-cumfreq[i75])/
    (cumfreq[i75+1]-cumfreq[i75]))
  lines(c(minx,end75),c(qnorm(.75,0,1),qnorm(.75,0,1)),lty=8)
  lines(c(end75,end75),c(-2.5,qnorm(.75,0,1)),lty=8)
  text(end75,-2.62,"Q3",cex=1.2)
}
