library(igraph)
library(spatstat)

dat <- read.csv("all_birds.csv")
dat <- na.omit(dat)
owin <- owin(xrange=c(min(dat$long),max(dat$long)), yrange=c(min(dat$lat),max(dat$lat)))
pp <- ppp(dat$long, dat$lat, window=owin, marks=dat$night)
pp <- unique.ppp(pp)
pp$marks <- as.factor(pp$marks)

pp1 <- pp[pp$marks=="1.1"|pp$marks=="1.2"]
pp2 <- pp[pp$marks=="2.1"|pp$marks=="2.2"]
pp3 <- pp[pp$marks=="3.1"|pp$marks=="3.2"]
nn1 <- nndist(pp1,  k=1, by=pp1$marks)
nn2 <- nndist(pp2,  k=1, by=pp2$marks)
nn3 <- nndist(pp3,  k=1, by=pp3$marks)
