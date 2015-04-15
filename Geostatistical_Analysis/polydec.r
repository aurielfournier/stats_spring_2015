# This function requires the following R-libraries.
library(tripack)
# =========================================================================================
# The "polydec" function takes a vectors of x- and y-coordinates and finds the influence
# polygons for each point, returning the areas of the polygons for use with the polygonal
# declustering global estimation method.  An outer edge is placed around the region at a
# distance roughly equal to half the minimum distance of a site to its neighbors, to handle
# edge effects.  To help place this boundary correctly, there is an optional "peels"
# argument that can be given.  The default is 1 meaning a single convex hull is taken.
# Generally, "peels" will not need to be any higher than 3 or 4.  This function calls the
# three functions in the "triangfuncs.r" file, so that this file must be compiled before
# executing "polydec".  Typing: "pd <- polydec(x,y,1)" will compute the influence polygons
# for all points in "x" and "y", returning the areas of these polygons in "pd".
# =========================================================================================
polydec <- function(x,y,peels=NA){
  if (missing(peels)) peels <- 1
  alldist <- dist(cbind(x,y))
  distmat <- matrix(0,nrow=attr(alldist,"Size"),ncol=attr(alldist,"Size"))
  distmat[lower.tri(distmat)] <- alldist
  distmat <- distmat + t(distmat)
  edgedist <- max(apply(distmat,1,function(a) a[order(a)[2]]))
  delaun <- triangles(tri.mesh(x,y))
  at <- t(delaun[,1:3])
  tx <- x; ty <- y
  ch <- rep(0,0)
  for (i in 1:peels){
    ch <- c(ch,chull(tx,ty))
    tx[ch] <- mean(tx)
    ty[ch] <- mean(ty)
  }
  atmat <- matrix(c(rep(at[1,],length(ch)),rep(at[2,],length(ch)),rep(at[3,],length(ch))),
    ncol=length(at[1,]),byrow=T)
  chmat <- matrix(rep(matrix(rep(ch,length(at[1,])),nrow=length(ch),byrow=F),3),ncol=length(at[1,]),
    byrow=F)
  equal <- atmat==chmat
  keep <- apply(equal,2,function(a) sum(a)>0)
  at <- at[,keep]
  pairs <- cbind(at[1:2,],at[c(1,3),],at[2:3,])
  extra.lab <- matrix(rep(pairs,6),nrow=2,byrow=F)
  pairdist <- sqrt((x[pairs[1,]]-x[pairs[2,]])^2+(y[pairs[1,]]-y[pairs[2,]])^2)
  pairs <- rbind(pairs,pairdist)
  end1.x <- (x[pairs[1,]]+x[pairs[2,]])/2 + (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end2.x <- (x[pairs[1,]]+x[pairs[2,]])/2 - (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end3.x <- x[pairs[1,]] + (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end4.x <- x[pairs[2,]] - (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end5.x <- x[pairs[1,]] - (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end6.x <- x[pairs[2,]] + (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end1.y <- (y[pairs[1,]]+y[pairs[2,]])/2 - (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end2.y <- (y[pairs[1,]]+y[pairs[2,]])/2 + (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end3.y <- y[pairs[1,]] - (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end4.y <- y[pairs[2,]] + (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end5.y <- y[pairs[1,]] + (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end6.y <- y[pairs[2,]] - (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  extra <- rbind(c(end1.x,end2.x,end3.x,end4.x,end5.x,end6.x),c(end1.y,end2.y,end3.y,end4.y,end5.y,end6.y))
  extra.out <- point.in.polygon(extra[1,],extra[2,],x[chull(x,y)],y[chull(x,y)])
  extra.x <- c(matrix(rep(extra[1,extra.out==F],2),nrow=2,byrow=F))
  extra.y <- c(matrix(rep(extra[2,extra.out==F],2),nrow=2,byrow=F))
  extra.lab <- c(t(extra.lab[,extra.out==F]))
  extra.nodes <- rbind(extra.x,extra.y,extra.lab)
  tr <- t(delaun[,1:3])
  an <- apply(tr,2,function(a,x,y) angles.tri(x[a],y[a]),x,y)

  maxan <- apply(an,2,max)
  an <- an[,maxan<180]
  tr <- tr[,maxan<180]
  acute <- apply(an,2,function(a) max(a)<90)
  tr.acute <- tr[,acute]
  mid.acute <- apply(tr.acute,2,function(a,x,y) mid.tri(x[a],y[a]),x,y)
  midseg1 <- apply(tr.acute,2,function(a,x,y) matrix(c(mean(x[a[1:2]]),
    mean(y[a[1:2]])),nrow=2,byrow=T),x,y)
  midseg2 <- apply(tr.acute,2,function(a,x,y) matrix(c(mean(x[a[c(1,3)]]),
    mean(y[a[c(1,3)]])),nrow=2,byrow=T),x,y)
  midseg3 <- apply(tr.acute,2,function(a,x,y) matrix(c(mean(x[a[2:3]]),
    mean(y[a[2:3]])),nrow=2,byrow=T),x,y)
  segs.ac <- cbind(rbind(x,y),matrix(rep(mid.acute,3),nrow=2),midseg1,
    midseg1,midseg2,midseg2,midseg3,midseg3)
  segs.ac <- rbind(segs.ac,c(1:length(x),c(t(tr.acute)),tr.acute[1,],tr.acute[2,],
    tr.acute[1,],tr.acute[3,],tr.acute[2,],tr.acute[3,]))

  an.ob <- an[,acute==F]
  tr.obtuse <- tr[,acute==F]
  mid.obtuse <- apply(tr.obtuse,2,function(a,x,y) mid.tri(x[a],y[a]),x,y)
  ext.nodes <- matrix(nrow=3,ncol=2*length(tr.obtuse[1,]))
  for (i in 1:length(tr.obtuse[1,])) ext.nodes[1:2,(2*i-1):(2*i)] <-
    edgenodes.tri(x[tr.obtuse[,i]],y[tr.obtuse[,i]],mid.obtuse[,i],an.ob[,i])
  ext.nodes <- cbind(ext.nodes,ext.nodes)
  obpoints0 <- apply(an.ob,2,function(a) order(a)[3])
  obpoints <- diag(apply(tr.obtuse,2,function(a,obpoints0) a[obpoints0],obpoints0))
  obpoints <- c(matrix(rep(obpoints,2),nrow=2,byrow=T))
  others1 <- rep(1,length(obpoints0))
  others1[obpoints0==1] <- 2
  others2 <- rep(3,length(obpoints0))
  others2[obpoints0==3] <- 2
  others1 <- diag(apply(tr.obtuse,2,function(a,others1) a[others1],others1))
  others2 <- diag(apply(tr.obtuse,2,function(a,others2) a[others2],others2))
  others <- c(matrix(c(others1,others2),nrow=2,byrow=T))
  ext.nodes[3,] <- c(obpoints,others)

  ob.in <- point.in.polygon(mid.obtuse[1,],mid.obtuse[2,],
    extra.x[chull(extra.x,extra.y)],extra.y[chull(extra.x,extra.y)])
  ob.x <- c(matrix(rep(mid.obtuse[1,ob.in==T],3),nrow=3,byrow=T))
  ob.y <- c(matrix(rep(mid.obtuse[2,ob.in==T],3),nrow=3,byrow=T))
  ob.lab <- c(tr.obtuse[,ob.in==T])
  ob.nodes <- rbind(ob.x,ob.y,ob.lab)

  mids.x <- c(matrix(rep((x[obpoints]+x[others])/2,2),nrow=2,byrow=T))
  mids.y <- c(matrix(rep((y[obpoints]+y[others])/2,2),nrow=2,byrow=T))
  mids.lab <- c(matrix(c(obpoints,others),nrow=2,byrow=T))
  mids.nodes <- rbind(mids.x,mids.y,mids.lab)

  allnodes <- cbind(segs.ac,ext.nodes,ob.nodes,mids.nodes,extra.nodes)
  dist.nodes <- sqrt((allnodes[1,]-x[allnodes[3,]])^2+(allnodes[2,]-y[allnodes[3,]])^2)
  min.dist <- apply(allnodes,2,function(a,x,y) min(sqrt((a[1]-x)^2 + (a[2]-y)^2)),x,y)
  minxydist <- min(dist(cbind(x,y)))
  omit.node <- abs(dist.nodes-min.dist)>(.00001*minxydist)
  allnodes <- allnodes[,omit.node==F]

  par(pty="s")
  xlims <- c(min(allnodes[1,]),max(allnodes[1,]))
  ylims <- c(min(allnodes[2,]),max(allnodes[2,]))
  plot(x,y,type="p",axes=T,xlab="",ylab="",xlim=xlims,ylim=ylims)
  polyarea <- rep(0,length(x))
  for (i in 1:length(x)){
    nodes <- allnodes[,allnodes[3,]==i]
    xord <- order(nodes[1,])
    polygon(nodes[1,chull(nodes[1,],nodes[2,])],nodes[2,chull(nodes[1,],nodes[2,])],
      density=0)
    polyarea[i] <- areapl(cbind(nodes[1,chull(nodes[1,],nodes[2,])],
      nodes[2,chull(nodes[1,],nodes[2,])]))
  }
  polyarea
}
