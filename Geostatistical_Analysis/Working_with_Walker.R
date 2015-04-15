# Creates a histogram of the Walker Lake V-values
# ===============================================
walk100 <- read.table("walk100.txt",header=T)   # Reads in 100-value Walker Lake data.
hbrk <- seq(-20,160,10)                              # Defines the histogram interval breaks.
hist(walk100$v+.1,breaks=hbrk,xlab="",ylab="",       # Creates a histogram of the V-values with no axis
     main="",xlim=c(-20,160),ylim=c(0,20),axes=F)       #   labels, axis limits set, and no axes plotted.
axis(1,at=seq(0,150,50),tick=F,pos=0,mgp=c(3,.5,0),  # Defines the x-axis(1), with axis labels at 0,50,
     cex.axis=1.5)                                      #   100, & 150 of character size 1.5.
axis(2,at=seq(0,20,5),pos=-20,tck=.02,mgp=c(3,0.5,0),# Defines the y-axis(2), with axis labels at 0,5,10,
     las=2,cex.axis=1.5)                                #   15, & 20 with inward tick marks (tck=.02).
mtext(side=1,line=2.0,"V (ppm)",cex=1.5)             # Puts an x-axis (side=1) label on the plot.
mtext(side=2,line=2.3,"Frequency (%)",cex=1.5)       # Puts a y-axis (side=2) label on the plot.

# Creates a histogram of the Walker Lake U-values
# ===============================================
hbrk <- seq(-10,70,5)
hist(walk100$u+.1,breaks=hbrk,xlab="",ylab="",main="",xlim=c(-10,70),ylim=c(0,20.1),axes=F)
axis(1,at=seq(0,80,20),tick=F,pos=0,mgp=c(3,.5,0),cex.axis=1.5)
axis(2,at=seq(0,20,5),pos=-10,tck=.02,mgp=c(3,0.5,0),las=2,cex.axis=1.5)
mtext(side=1,line=2.0,"U (ppm)",cex=1.5)
mtext(side=2,line=2.3,"Frequency (%)",cex=1.5)
lines(c(15,17.5),c(20.9,20)); lines(c(17.5,20),c(20,20.9)); text(17.5,19.4,"33",cex=1.5)

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
ppnorm(walk100$v,seq(10,140,10),c(0,50,150),"V (ppm)") # Creates a normal PP-plot of the V-data.
ppnorm(walk100$u,seq(5,40,5),c(0,5,30,40),"U (ppm)")   # Creates a normal PP-plot of the U-data.

# Creates a text plot of the V- and U-values
# ==========================================
plot(walk100$x,walk100$y,pch="+",xlab="",ylab="",  # Plots the (x,y) locations with a "+" sign, no
     xlim=c(.7,10.3),ylim=c(.3,10.7),axes=F)          #   axis labels, specified axis limits, & no axes.
box()                                              # Places a box around the plot.
for (i in 1:100){                                  # Begins a for loop.
  text(walk100$x[i],walk100$y[i]+.30,              # Iteratively writes the V-values over the "+"
       paste(walk100$v)[i],cex=0.9)                   #   signs of character size 0.9 (default=1).
  text(walk100$x[i],walk100$y[i]-.25,              # Iteratively writes the U-values under the "+"
       paste(walk100$u)[i],cex=0.9)                   #   signs of character size 0.9.
}                                                  # Ends the for loop.

# Creates a quantile-quantile (Q-Q) plot of the V and U-values
# ============================================================
uq <- quantile(walk100$u,probs=seq(.05,.95,.05))   # Computes the .05, .10, ..., .95 quantiles of U.
vq <- quantile(walk100$v,probs=seq(.05,.95,.05))   # Computes the .05, .10, ..., .95 quantiles of V.
plot(vq,uq,pch="+",xlab="",ylab="",axes=F,         # Plots the U-quantiles vs. V-quantiles, with no
     xlim=c(0,150),ylim=c(0,60),cex=1.5)              #   axis labels, no axes, and axis limits set.
axis(1,at=seq(0,150,50),tck=.02,pos=0,             # Puts an x-axis on the plot, at y=0 (pos=0), with
     mgp=c(3,.5,0),cex.axis=1.5)                      #   inward ticks (tck=.02) at 0,50,100,150 (at).
axis(2,at=seq(0,60,20),tck=.02,pos=0,              # Puts a y-axis on the plot, at x=0 (pos=0), with
     mgp=c(3,0.5,0),las=2,cex.axis=1.5)               #   inward ticks (tck=.02) at 0,20,40,60 (at).
mtext(side=1,line=2.3,"V Quantiles (ppm)",cex=1.5) # Puts an x-axis (side=1) label on the plot.
mtext(side=2,line=2.0,"U Quantiles (ppm)",cex=1.5) # Puts a y-axis (side=2) label on the plot.
lines(c(0,60),c(0,60),lty=8)                       # Plots dashed (lty=8) line from (0,0) to (60,60).
text(50,40,"U=V",cex=1.5)                          # Puts text "U=V" on the plot at (x,y)=(50,40).
lines(c(vq[5],88),c(uq[5],4))                      # The next 6 lines put text labels on the plot
text(115,4,"Lower Quartile (81.3,14.0)")           #   with lines drawn to them.
lines(c(vq[10],96),c(uq[10],9))
text(115,8,"Median (100.5,18.0)")
lines(c(vq[15],110),c(uq[15],13))
text(128,12,"Upper Quartile (116.8,25.0)")

# Creates a scatterplot of the V-values vs. U-values with the regression line
# ===========================================================================
plot(walk100$v,walk100$u,pch="+",xlab="",ylab="",            # Plots U vs. V with the "+" sign, no axis
     axes=F,xlim=c(0,150),ylim=c(0,60),cex=1.5)                 #   labels, no axes, and axis limits as set.
axis(1,at=seq(0,150,50),tck=.02,pos=0,mgp=c(3,.5,0),         # Puts an x-axis on the plot.
     cex.axis=1.5)
axis(2,at=seq(0,60,20),tck=.02,pos=0,mgp=c(3,0.5,0),las=2,   # Puts a y-axis on the plot.
     cex.axis=1.5)
mtext(side=1,line=2.3,"V (ppm)",cex=1.5)                     # Puts an x-axis label on the plot.
mtext(side=2,line=2.0,"U (ppm)",cex=1.5)                     # Puts a y-axis label on the plot.
uv.reg <- lsfit(walk100$v,walk100$u)                         # Performs least squares regression of U on V.
abline(uv.reg,lwd=2)                                         # Overlays the regression line on the plot.

# Creates a scatterplot of the V-values vs. U-values with a loess fit plotted
# ===========================================================================
plot(walk100$v,walk100$u,pch="+",xlab="",ylab="",            # Plots U vs. V with the "+" sign, no axis
     axes=F,xlim=c(0,150),ylim=c(0,60),cex=1.5)                 #   labels, no axes, and axis limits as set.
axis(1,at=seq(0,150,50),tck=.02,pos=0,mgp=c(3,.5,0),         # Puts an x-axis on the plot.
     cex.axis=1.5)
axis(2,at=seq(0,60,20),tck=.02,pos=0,mgp=c(3,0.5,0),las=2,   # Puts a y-axis on the plot.
     cex.axis=1.5)
mtext(side=1,line=2.3,"V (ppm)",cex=1.5)                     # Puts an x-axis label on the plot.
mtext(side=2,line=2.0,"U (ppm)",cex=1.5)                     # Puts a y-axis label on the plot.
lines(lowess(walk100$v,walk100$u,f=1/3),lwd=2)               # Puts a loess fit through the data where
#   1/3 the values are used for local fits.
# Creates text plots of the low and high V-values
# ===============================================
par(mfrow=c(1,1))                                  # Creates a 1x1 graphics window.
plot(walk100$x,walk100$y,pch="+",xlab="",ylab="",  # Plots the (x,y) locations with a "+" sign, no
     xlim=c(.5,10.5),ylim=c(.5,11),axes=F,cex=1.3)    #   axis labels, specified axis limits, & no axes.
box()                                              # Places a box around the plot.
text(0.5,11,"(a)",cex=1.5)                         # Puts text on the plot at (.5,11).
for (i in 1:100){                                  # Begins a for loop.
  text(walk100$x[i],walk100$y[i]+.4,               # Iteratively writes the V-values over the "+"
       paste(walk100$v)[i],cex=1.5)                   #   signs of character size 0.9 (default=1).
  if (walk100$v[i]<=70) polygon(walk100$x[i]+      # Puts a box around all V-values which are 70
                                  c(-.29,-.29,.29,.29),walk100$y[i]+             #   or less.
                                  c(-.2,.7,.7,-.2),density=0)
}                                                  # Ends the for loop.
plot(walk100$x,walk100$y,pch="+",xlab="",ylab="",  # Plots the (x,y) locations with a "+" sign, no
     xlim=c(.5,10.5),ylim=c(.5,11),axes=F,cex=1.3)    #   axis labels, specified axis limits, & no axes.
box()                                              # Places a box around the plot.
text(0.5,11,"(b)",cex=1.5)                         # Puts text on the plot at (.5,11).
for (i in 1:100){                                  # Begins a for loop.
  text(walk100$x[i],walk100$y[i]+.4,               # Iteratively writes the V-values over the "+"
       paste(walk100$v)[i],cex=1.5)                   #   signs of character size 0.9 (default=1).
  if (walk100$v[i]>=128) polygon(walk100$x[i]+     # Puts a box around all V-values which are 128
                                   c(-.38,-.38,.38,.38),walk100$y[i]+             #   or more.
                                   c(-.2,.7,.7,-.2),density=0)
}                                                  # Ends the for loop.

# Creates a contour plot of the V-values (n=100)
# ==============================================
library(akima)                                   # Loads "akima" library for interpolation.
int.v <- interp(walk100$x,walk100$y,             # Interpolates the V-values over the 10x10
                walk100$v)                                     #   region.
contour(int.v,levels=seq(0,140,5),xlab="",       # Creates a contour plot with the interpolated
        ylab="",xlim=c(1.7,9.7),ylim=c(1.2,9.7),       #   v-values at levels 0,5,10, ..., 140.
        labcex=1,axes=F)
box()                                            # Puts a box around the plot.

# Creates a symbol map of the 100 V-data
# ======================================
sym <- trunc(walk100$v/15)                           # Computes the symbols for each V-value.
plot(walk100$x,walk100$y,pch=" ",xlab="",ylab="",    # Makes a blank plot of the V-locations, with
     xlim=c(.5,15.5),ylim=c(.5,11),axes=F)              #   no axes or axis labels, and axis limits.
polygon(c(.5,.5,10.5,10.5),c(.35,10.7,10.7,.35),     # Places a box around the plot without filling
        density=0)                                         #   it in (density=0).
for (i in 1:100){                                    # Begins a for loop.
  text(walk100$x[i],walk100$y[i],                    # Iteratively writes the V-symbols at their
       paste(sym)[i],cex=1.5)                           #   locations of size 1.5 (default=1).
}                                                    # Ends the for loop.
labs <- c("0 = 0-14 ppm","1 = 15-29","2 = 30-44",    # Assigns the labels for the symbol legend.
          "3 = 45-59","4 = 60-74","5 = 75-89","6 = 90-104",
          "7 = 105-119","8 = 120-134","9 = 135-149")
legend(10.2,9,labs,bty="n",cex=1.5)                  # Puts a legend on the plot with no box (bty).

# Creates a greyscale plot of the V-values (n=100)
# ================================================
image(1:10,1:10,t(matrix(walk100$v,nrow=10)),  # Creates a greyscale image plot on the 10x10 grid,
      xlim=c(1,14),xlab="",ylab="",axes=F,         #   with extra space at the right (xlim goes to 14),
      col=rev(heat.colors(24)))                    #   and reversed topographic colors 24-1 .
# REMEMBER TO RUN THE IMAGE.LEGEND FUNCTION HERE!!!!!!!
image.legend(11,9,zlim=range(walk100$v),       # Puts a legend on the map with upper left corner
             col=rev(heat.colors(24)))                    #   at (11,9) according to the 24 colors.

# Creates a sequence of indicator maps of the V-values (n=100)
# ============================================================
par(mfrow=c(3,3),pty="m",mar=c(1,2,1,2))      # Creates a 3x3 plotting grid with reduced margins.
cutoffs <- seq(15,135,15)                     # Sets the indicator plot cutoff levels.
plotlabs <- c("(a)","(b)","(c)","(d)","(e)",  # Creates a vector of labels for the 9 plots.
              "(f)","(g)","(h)","(i)")
for (i in 1:length(cutoffs)){                 # Begin for loop for the 9 plots.
  labs <- c(paste("V <",cutoffs[i],"ppm"),    # Defines the 2 labels for the legend.
            paste("V >=",cutoffs[i],"ppm"))
  image(1:10,1:10,t(matrix(walk100$v>=cutoffs[i], # Creates an indicator plot according to the
                           nrow=10)),ylim=c(-2.2,14),axes=F,xlab="",     #   cutoff specified, with no axes, and 2 colors
        ylab="",col=c(0,1))                           #   (0=white, 1=black).
  legend(0,14,labs,density=c(0,-1),bty="n")       # Puts a (white,black)=(0,-1) legend on the plot.
  polygon(c(.5,.5,10.5,10.5),c(.5,10.5,10.5,.5),  # Draws a box around the plot, and does not fill
          density=0)                                    #   it in (density=0).
  text(5.5,-1,plotlabs[i],cex=1.5)                # Puts the text label on the plot at (5.5,-1).
}                                             # End of for loop.

# Creates text plots of V-values to show moving window calculations
# =================================================================
par(mfrow=c(1,1))
plot(walk100$x,walk100$y,pch="+",xlab="",ylab="",  # Plots the (x,y) locations with a "+" sign, no
     xlim=c(.5,10.5),ylim=c(.5,11),axes=F)            #   axis labels, specified axis limits, & no axes.
box()                                              # Places a box around the plot.
polygon(c(.5,.5,4.5,4.5),c(6.7,10.7,10.7,6.7),     # Puts a 4x4 window in the upper LH corner of
        density=0,lty=8)                                 #   the text plot using a dashed line (lty=8).
for (i in 1:100){                                  # Begins a for loop.
  text(walk100$x[i],walk100$y[i]+.3,               # Iteratively writes the V-values over the "+"
       paste(walk100$v)[i],cex=1.5)                   #   signs of character size 0.9 (default=1).
}                                                  # Ends the for loop.
text(0.4,11.1,"(a)",cex=1.5)                       # Puts the text "(a)" on the plot.
plot(walk100$x,walk100$y,pch="+",xlab="",ylab="",  # Plots the (x,y) locations with a "+" sign, no
     xlim=c(.5,10.5),ylim=c(.5,11),axes=F)            #   axis labels, specified axis limits, & no axes.
box()                                              # Places a box around the plot.
polygon(c(2.5,2.5,6.5,6.5),c(6.7,10.7,10.7,6.7),   # Puts a 4x4 window 2 units right of the upper
        density=0,lty=8)                                 #   LH corner using a dashed line (lty=8).
for (i in 1:100){                                  # Begins a for loop.
  text(walk100$x[i],walk100$y[i]+.3,               # Iteratively writes the V-values over the "+"
       paste(walk100$v)[i],cex=1.5)                   #   signs of character size 0.9 (default=1).
}                                                  # Ends the for loop.
text(0.4,11.1,"(b)",cex=1.5)                       # Puts the text "(b)" on the plot.

# Calculates moving window means and SDs for the Walker Lake data
# ===============================================================
# REMEMBER TO RUN THE MOVEWIN FUNCTION HERE!!!!!!!
move100 <- movewin(walk100$x,walk100$y,walk100$v,  # Computes moving windows statistics using 4x4
                   4,4,2)                                           #   windows with 2 units of overlap.

# Creates a text plot of the moving windows means and SDs
# =======================================================
dox <- min(move100$x)-.5                           # Defines the lower x-axis limit.
doy <- min(move100$y)-.5                           # Defines the lower y-axis limit.
upx <- max(move100$x)+.5                           # Defines the upper x-axis limit.
upy <- max(move100$y)+.5                           # Defines the upper y-axis limit.
plot(move100$x,move100$y,pch="+",xlab="",ylab="",  # Plots the (x,y) locations with a "+" sign, no
     xlim=c(dox,upx),ylim=c(doy,upy),axes=F,cex=1.3)  #   axis labels, specified axis limits, & no axes.
box()                                              # Places a box around the plot.
for (i in 1:100){                                  # Begins a for loop.
  text(move100$x[i],move100$y[i]+.25,              # Iteratively writes the V-values over the "+"
       format(round(move100$means,1))[i],cex=1.3)     #   signs of character size 1.3 (default=1).
  text(move100$x[i],move100$y[i]-.25,              # Iteratively writes the U-values under the "+"
       format(round(move100$sdevs,1))[i],cex=1.3)     #   signs of character size 1.3.
}                                                  # Ends the for loop.

# Creates a plot of the moving window standard deviations vs. means
# =================================================================
par(mar=c(5, 4, 4, 2) + 0.1)                       # Resets margin widths to default settings.
plot(move100$means,move100$sdevs,pch="+",xlab="",    # Plots the moving window SD's vs. means, with
     ylab="",axes=F,xlim=c(80,110),ylim=c(0,50),cex=1.5)#   no axes, and plotting symbol "+".
axis(1,at=seq(80,110,10),tck=.02,pos=0,            # Puts an x-axis on the plot, at y=0 (pos=0), with
     mgp=c(3,.5,0),cex.axis=1.5)                      #   inward ticks (tck=.02) at 0,50,100,150 (at).
axis(2,at=seq(0,50,10),tck=.02,pos=80,             # Puts a y-axis on the plot, at x=0 (pos=0), with
     mgp=c(3,0.5,0),las=2,cex.axis=1.5)               #   inward ticks (tck=.02) at 0,20,40,60 (at).
mtext(side=1,line=2.0,"Mean",cex=1.8)              # Puts an x-axis (side=1) label on the plot.
mtext(side=2,line=2.0,"Standard Deviation",cex=1.8)# Puts a y-axis (side=2) label on the plot.

# Creates h-scatterplots at orientations (0,1)--(0,4) for V-values
# ================================================================
par(mfrow=c(2,2),pty="m",mar=c(3.6,3.6,3.6,2.1))# Creates a 2x2 plotting grid with reduced margins.
V <- walk100$v; x <- walk100$x; y <- walk100$y  # Defines V,x,y to have shorter names.
hscatter(x,y,V,V,h=c(0,1))                      # Creates an h=(0,1) scatterplot of V(t+h) vs. V(t).
hscatter(x,y,V,V,h=c(0,2))                      # Creates an h=(0,2) scatterplot of V(t+h) vs. V(t).
hscatter(x,y,V,V,h=c(0,3))                      # Creates an h=(0,3) scatterplot of V(t+h) vs. V(t).
hscatter(x,y,V,V,h=c(0,4))                      # Creates an h=(0,4) scatterplot of V(t+h) vs. V(t).

# _____________________________________________________________________________________________________________
# Plots the covariogram, correlogram, and semivariogram for the
# Walker Lake (n=100) V-data.  The values were determined using
# the hscatter function found on the course webpage.
# =============================================================
h <- c(1,2,3,4,5,6)                           # Defines the six h-values used.
ch <- c(448.8,341,323.8,291.5,247.5,221.6)    # Defines the vector of covariances.
ph <- c(.742,.590,.560,.478,.365,.285)        # Defines the vector of correlogram values.
gh <- c(156.4,239.6,260.7,326.5,440.9,571.4)  # Defines the vector of semivariogram values.
plot(h,ch,type="l",xlab="|h|",ylab="C(h)",    # Plots the covariance function, with x-axis
     cex.lab=1.6,cex.axis=1.5,cex.main=1.8,      #   label, y-axis label, and title, printed
     main="Covariogram for V",mgp=c(2.7,1,0))    #   according to the sizes in "cex".
plot(h,ph,type="l",xlab="|h|",ylab="p(h)",    # Plots the covariance function, with x-axis
     cex.lab=1.6,cex.axis=1.5,cex.main=1.8,      #   label, y-axis label, and title, printed
     main="Correlogram for V",mgp=c(2.7,1,0))    #   according to the sizes in "cex".
plot(h,gh,type="l",xlab="|h|",ylab="Gamma(h)",# Plots the covariance function, with x-axis
     cex.lab=1.6,cex.axis=1.5,cex.main=1.8,      #   label, y-axis label, and title, printed
     main="Semivariogram for V",mgp=c(2.7,1,0))  #   according to the sizes in "cex".



# _____________________________________________________________________________________________________________
# Reads in Walker Lake Sample Data Set (n=470)
# ============================================
walk470 <- read.table("walk470.txt",header=T)   # Reads in the data.
par(mfrow=c(1,1))                                  # Creates a 1x1 graphics window.

plot(walk470$x,walk470$y,xlab="x",ylab="y",          # Plots the 470 (x,y)-locations of the sample
     cex.lab=1.6,cex.axis=1.5,cex.main=1.8,             #   Walker Lake data.
     main="Walker Lake Sample Locations\n(U-locations are closed circles)")
points(walk470$x[!is.na(walk470$u)],                 # Overlays the 295 U-locations as closed
       walk470$y[!is.na(walk470$u)],pch=16)               #   circles.

# Creates a histogram of the Walker Lake Sample V-values
# ======================================================
hbrk <- seq(-100,1550,50)                            # Defines the histogram interval breaks.
hist(walk470$v+.1,breaks=hbrk,xlab="V (ppm)",ylab=   # Creates a histogram of the V-values with axis
       "Frequency (%)",xlim=c(-100,1550),ylim=c(0,70.5),  #   labels, axis limits set, and no axes plotted.
     main="",cex.lab=1.6,axes=F)
axis(1,at=seq(0,1500,500),tick=F,pos=0,              # Defines the x-axis(1), with axis labels at 0,500,
     mgp=c(3,.5,0),cex.axis=1.5)                        #   1000, & 1500 of character size 1.5.
axis(2,at=seq(0,15,5)*4.7,pos=-100,mgp=c(3,0.5,0),   # Defines the y-axis(2), with axis labels at 0,5,10,
     tck=.02,labels=seq(0,15,5),las=2,cex.axis=1.5)     #   & 15% with inward tick marks (tck=.02).

# Creates a histogram of the Walker Lake Sample U-values
# ======================================================
hbrk <- seq(-100,1550,50)                            # Defines the histogram interval breaks.
u <- walk470$u[walk470$u<1550]                       # Eliminates U-values at 1550 or above.
hist(u+.1,breaks=hbrk,xlab="U (ppm)",ylab=           # Creates a histogram of the U-values with axis
       "Frequency (%)",xlim=c(-100,1550),ylim=c(0,55),    #   labels, axis limits set, and no axes plotted.
     main="",cex.lab=1.6,axes=F)
axis(1,at=seq(0,1500,500),tick=F,pos=0,              # Defines the x-axis(1), with axis labels at 0,500,
     mgp=c(3,.5,0),cex.axis=1.5)                        #   1000, & 1500 of character size 1.5.
axis(2,at=seq(0,20,5)*2.75,pos=-100,mgp=c(3,0.5,0),  # Defines the y-axis(2), with axis labels at 0,5,10,
     tck=.02,labels=seq(0,20,5),las=2,cex.axis=1.5)     #   15 & 20% with inward tick marks (tck=.02).

# Creates a scatterplot of the Sample U-values vs. V-values
# =========================================================
type <- walk470$t
plot(walk470$v,walk470$u,pch=" ",xlab="V (ppm)",ylab="U (ppm)", # Plots U vs. V with the "+" sign, no axis
     axes=F,xlim=c(0,1500),ylim=c(0,2000),cex.lab=1.6)             #   labels, no axes, and axis limits as set.
points(walk470$v[type==1],walk470$u[type==1],pch="o",cex=1.5)   # Plots the type 1 points.
points(walk470$v[type==2],walk470$u[type==2],pch="+",cex=1.5)   # Plots the type 2 points.
axis(1,at=seq(0,1500,500),tck=.02,pos=0,mgp=c(3,.5,0),          # Puts an x-axis on the plot, with inward
     cex.axis=1.3)                                                 #   tick marks, at position 0.
axis(2,at=seq(0,2000,500),tck=.02,pos=0,mgp=c(3,.5,0),las=2,cex.axis=1.3) # Puts a y-axis on the plot.

# Creates a plot of moving window means vs. SDs for V data and U data
# ===================================================================
v.move <- movewin(walk470$x,walk470$y,walk470$v,60,60,40)       # Calculates moving window means & SDs
plot(v.move$means,v.move$sdevs,xlim=c(0,1000),ylim=c(0,500),    # Plots the moving window SDs vs. means
     pch=" ",xlab="Mean",ylab="Standard Deviation",axes=F,         #   for the V-data with axis labels and
     cex.lab=1.6,cex=1.5)                                          #   a blank plotting character.
points(v.move$means[v.move$numvals<20],                         # Plots (SD,mean) pairs for windows with
       v.move$sdevs[v.move$numvals<20],pch="o",cex=1.5)              #   less than 20 points.
points(v.move$means[v.move$numvals>=20],                        # Plots (SD,mean) pairs for windows with
       v.move$sdevs[v.move$numvals>=20],pch="+",cex=1.5)             #   20 points or more.
axis(1,at=seq(0,1000,200),tck=.02,pos=0,mgp=c(3,.5,0),          # Puts an x-axis on the plot.
     cex.axis=1.5)
axis(2,at=seq(0,500,100),tck=.02,pos=0,mgp=c(3,.5,0),las=2,     # Puts a y-axis on the plot.
     cex.axis=1.5)
text(100,510,"(a)",cex=1.5)                                     # Puts text on plot at (100,510).

u.move <- movewin(walk470$x[196:470],walk470$y[196:470],        # Calculates moving window means & SDs
                  walk470$u[196:470],60,60,40)                                  #   for U-data without missing values.
plot(u.move$means,u.move$sdevs,xlim=c(0,1500),ylim=c(0,1500),   # Plots the moving window SDs vs. means
     pch=" ",xlab="Mean",ylab="Standard Deviation",axes=F,         #   for the U-data with axis labels and
     cex.lab=1.6,cex=1.5)                                          #   a blank plotting character.
points(u.move$means[u.move$numvals<20],                         # Plots (SD,mean) pairs for windows with
       u.move$sdevs[u.move$numvals<20],pch="o",cex=1.5)              #   less than 20 points.
points(u.move$means[u.move$numvals>=20],                        # Plots (SD,mean) pairs for windows with
       u.move$sdevs[u.move$numvals>=20],pch="+",cex=1.5)             #   20 points or more.
axis(1,at=seq(0,1500,500),tck=.02,pos=0,mgp=c(3,.5,0),          # Puts an x-axis on the plot.
     cex.axis=1.5) 
axis(2,at=seq(0,1500,500),tck=.02,pos=0,mgp=c(3,1.0,0),         # Puts a y-axis on the plot.
     las=2,cex.axis=1.5)
text(120,1520,"(b)",cex=1.5)                                   # Puts text on plot at (-170,1650).

# Creates a plot of moving window means vs. SDs for log-transformed V data and U data
# ===================================================================================
logv <- log(walk470$v+1)                                        # Log transform of V-values.
v.move <- movewin(walk470$x,walk470$y,logv,60,60,40)            # Calculates moving window means & SDs
plot(v.move$means,v.move$sdevs,xlim=c(0,7),ylim=c(0,3),pch=" ", # Plots the moving window SDs vs.
     xlab="Mean (Log(V+1))",ylab="Standard Deviation",axes=F,      #   means for the log V-data with
     main="Means vs. SDs for Log-Transformed V-Data",cex.lab=1.6,  #   title and axis labels.
     cex.main=1.6,cex=1.5)
points(v.move$means[v.move$numvals<20],                         # Plots (SD,mean) pairs for windows with
       v.move$sdevs[v.move$numvals<20],pch="o",cex=1.5)              #   less than 20 points.
points(v.move$means[v.move$numvals>=20],                        # Plots (SD,mean) pairs for windows with
       v.move$sdevs[v.move$numvals>=20],pch="+",cex=1.5)             #   20 points or more.
axis(1,at=seq(0,7,1),tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=1.5)  # Puts an x-axis on the plot.
axis(2,at=seq(0,3,1),tck=.02,pos=0,mgp=c(3,.5,0),las=2,         # Puts a y-axis on the plot.
     cex.axis=1.5)

logu <- log(walk470$u[196:470]+1)                               # Log transform of U-values.
u.move <- movewin(walk470$x[196:470],walk470$y[196:470],        # Calculates moving window means & SDs
                  logu,60,60,40)                                                #   for U-data without missing values.
plot(u.move$means,u.move$sdevs,xlim=c(0,7),ylim=c(0,5),pch=" ", # Plots the moving window SDs vs.
     xlab="Mean (Log(U+1))",ylab="Standard Deviation",,axes=F,     #   means for the log U-data with
     main="Means vs. SDs for Log-Transformed U-Data",cex.lab=1.6,  #   title and axis labels.
     cex.main=1.6,cex=1.5)
points(u.move$means[u.move$numvals<20],                         # Plots (SD,mean) pairs for windows with
       u.move$sdevs[u.move$numvals<20],pch="o",cex=1.5)              #   less than 20 points.
points(u.move$means[u.move$numvals>=20],                        # Plots (SD,mean) pairs for windows with
       u.move$sdevs[u.move$numvals>=20],pch="+",cex=1.5)             #   20 points or more.
axis(1,at=seq(0,7,1),tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=1.5)  # Puts an x-axis on the plot.
axis(2,at=seq(0,5,1),tck=.02,pos=0,mgp=c(3,.5,0),las=2,         # Puts a y-axis on the plot.
     cex.axis=1.5)

# Creates a plot of moving window means vs. SDs for square root-transformed V data and U data
# ===========================================================================================
sqrtv <-sqrt(walk470$v)                                         # Square root transform of V-values.
v.move <- movewin(walk470$x,walk470$y,sqrtv,60,60,40)           # Calculates moving window means & SDs
plot(v.move$means,v.move$sdevs,xlim=c(0,26),ylim=c(0,13),       # Plots the moving window SDs vs.
     xlab="Mean (Square Root(V))",ylab="Standard Deviation",       #   means for the square root V-data
     main="Means vs. SDs for Square Root-Transformed V-Data",
     pch=" ",axes=F,cex.lab=1.6,cex.main=1.6,cex=1.5)
points(v.move$means[v.move$numvals<20],                         # Plots (SD,mean) pairs for windows with
       v.move$sdevs[v.move$numvals<20],pch="o",cex=1.5)              #   less than 20 points.
points(v.move$means[v.move$numvals>=20],                        # Plots (SD,mean) pairs for windows with
       v.move$sdevs[v.move$numvals>=20],pch="+",cex=1.5)             #   20 points or more.
axis(1,at=seq(0,25,5),tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=1.5) # Puts an x-axis on the plot.
axis(2,at=seq(0,12,3),tck=.02,pos=0,mgp=c(3,.5,0),las=2,        # Puts a y-axis on the plot.
     cex.axis=1.5)

sqrtu <- sqrt(walk470$u[196:470])                               # Square root transform of U-values.
u.move <- movewin(walk470$x[196:470],walk470$y[196:470],        # Calculates moving window means & SDs
                  sqrtu,60,60,40)                                               #   for sqrt U-data w/o missing values.
plot(u.move$means,u.move$sdevs,xlim=c(0,32),ylim=c(0,18),       # Plots the moving window SDs vs.
     xlab="Mean (Square Root(U))",ylab="Standard Deviation",       #   means for the square root U-data
     main="Means vs. SDs for Square Root-Transformed U-Data",
     pch=" ",axes=F,cex.lab=1.6,cex.main=1.6,cex=1.5)
points(u.move$means[u.move$numvals<20],                         # Plots (SD,mean) pairs for windows with
       u.move$sdevs[u.move$numvals<20],pch="o",cex=1.5)              #   less than 20 points.
points(u.move$means[u.move$numvals>=20],                        # Plots (SD,mean) pairs for windows with
       u.move$sdevs[u.move$numvals>=20],pch="+",cex=1.5)             #   20 points or more.
axis(1,at=seq(0,32,4),tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=1.5) # Puts an x-axis on the plot.
axis(2,at=seq(0,18,3),tck=.02,pos=0,mgp=c(3,.5,0),las=2,        # Puts a y-axis on the plot.
     cex.axis=1.5)

library(tripack)                # Needed for "tri.mesh" triangulation
library(sp)                     # Needed for "point.in.polygon" function
library(splancs)                # Needed for "areapl" function
source("triangfuncs.r") # Loads triangulation functions for use
#   with polydec.
source("allfunctions.r")# Loads R-script with all of my functions

# Demo with small (n=10) data from page 49 of class notes
# =======================================================
x <- c(0,4,5,8,9,11.5,11.9,12,12.7,13)
y <- c(1.3,3.5,.4,3,0,2.2,3.9,3,2.4,3.3)*4
polydec(x,y,peels=1)
text(0,4.4,"V1"); text(4,14.8,"V2"); text(5,0.9,"V3")
text(8,12.7,"V4"); text(9,-.6,"V5"); text(11.5,8.2,"V6")
text(12.0,16.3,"V7"); text(11.2,12.0,"V8"); text(13.6,9.6,"V9")
text(14.0,13.2,"V10")

walk470 <- read.table("walk470.txt",header=T)

# Performs polygonal declustering and returns the site weights
# ============================================================
walk.pd <- polydec(walk470$x,walk470$y,peels=5)
title("Polygonal Declustering",cex.main=1.6)

# Performs cell declustering for a given window size
# ==================================================
xval <- 200/seq(5,50,.5); yval <- 300/seq(5,50,.5)
globmean <- matrix(nrow=length(xval),ncol=length(yval))
for (i in 1:length(xval)){
  for (j in 1:length(yval)){
    globmean[i,j] <- celldec(walk470$x,walk470$y,walk470$v,i,j)
  }
}
library(akima)
V.int <- interp(rep(xval,length(yval)),rep(yval,each=length(xval)),
                c(globmean),xo=seq(5,50,0.5),yo=seq(5,50,0.5))
contour(V.int,xlim=c(5,40))


# Plots the cross h-scatterplots for the Walker Lake (n=100) V&U-data, and
# plots the cross-covariogram, cross-correlogram, and cross-semivariogram
# functions for the U,V data.  The hscatter function was used to generate
# the cross h-scatterplots.
# ========================================================================
walk100 <- read.table("Data/walk100.txt",header=T) # Reads in walk100 Walker Lake data.
x <- walk100$x                         # x is set to the x-values of "walk100".
y <- 11 - walk100$y                    # y is set to 11 - the y-values of "walk100".
v <- walk100$v                         # v is set to the V-values of "walk100".
u <- walk100$u                         # u is set to the U-values of "walk100".

# Creates the 5 cross-h scatterplots on page 39 of the class notes
# ================================================================
par(mfrow=c(3,2))                      # Sets up a 3x2 graphics window.
hscatter(x,y,v,u,c(0,0))               # Produces a scatterplot of u vs. v.
hscatter(x,y,v,u,c(0,1))               # Produces a cross h=(0,1)-scatterplot of u vs. v.
hscatter(x,y,v,u,c(0,2))               # Produces a cross h=(0,2)-scatterplot of u vs. v.
hscatter(x,y,v,u,c(0,3))               # Produces a cross h=(0,3)-scatterplot of u vs. v.
hscatter(x,y,v,u,c(1,0))               # Produces a cross h=(1,0)-scatterplot of u vs. v.

# Plots the 3 cross-functions on page 42 of the class notes
# =========================================================
h <- c(0,1,2,3,4,5,6)                              # Sets a vector of h-values.
out <- matrix(nrow=length(h),ncol=3)               # Defines an hx3 blank matrix.
for (i in 0:6) out[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,v,u,c(0,i)))                        #   calculates cross-functions for each.
cch <- out[,1]; cph <- out[,2]; cgh <- out[,3]     # Defines the vectors of cross-covariances,
#   cross-corr's, & cross-semivariograms.
par(mfrow=c(2,2))                                  # Sets up a 2x2 graphics window.
plot(h,cch,type="n",axes=F,xlab="",ylab="")        # Creates a completely blank plot.
plot(h,cch,type="l",xlab="|h|",ylab="Cross - C(h)",# Plots the cross-correlations vs. h
     cex.lab=1.5,cex.axis=1.3)                        #   with sizes controlled by "cex".
title("Cross-Covariogram for U,V",cex=1)           # Puts a title on the plot.
plot(h,cph,type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     cex.lab=1.5,cex.axis=1.3)                        #   with sizes controlled by "cex".
title("Cross-Correlogram for U,V",cex=1)           # Puts a title on the plot.
plot(h[-1],cgh[-1],type="l",xlab="|h|",cex.lab=1.5,# Plots the cross-semivariogram vs. h
     ylab="Cross - Gamma(h)",cex.axis=1.3)            #   with sizes controlled by "cex".
title("Cross-Semivariogram for U,V",cex=1)         # Puts a title on the plot.


