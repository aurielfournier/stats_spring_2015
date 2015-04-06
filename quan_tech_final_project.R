library(ggplot2)
library(rgdal)
library(gridExtra)

dat <- read.csv("all_birds.csv")
dat <- na.omit(dat)

utm <- as.data.frame(project(cbind(dat$long, dat$lat), "+proj=utm +zone=15 ellps=WGS84"))
colnames(utm) <- c("utm_w","utm_n")

dat <- cbind(dat, utm)
dat4 <- dat[dat$year==2014,]


jdate4 <- unique(dat4$jdate)


list4 <- list()


## so this takes each day adn makes a list with each unique date being a level in the lsit
# removes days in which only one observer saw birds
for(i in 1:length(jdate4)){
  dat <- dat4[dat4$jdate==jdate4[i],]
  if (length(unique(dat$obs))>1){
  list4[[i]] <- dat}
}



## this pairs up eahc point with the nearest point by the other observer and measures the difference

dist4 <- list()
newdf <- list()

for (i in c(1:14,16:21,24:40)) {
  df <- as.data.frame(list4[[i]])
  if (nrow(df)>1){
  a <- df[df$obs=="N",]
  b <- df[df$obs=="A",]
            for(j in 1:nrow(a)){
            c <- b[which.min(sqrt(abs(a[j,"utm_n"]-b[,"utm_n"])^2 + abs(a[j,"utm_w"]-b[,"utm_w"])^2)),]
            d <- a[j,]
            c$dist <- min(sqrt(abs(a[j,"utm_n"]-b[,"utm_n"])^2 + abs(a[j,"utm_w"]-b[,"utm_w"])^2))
            dist4[[j]] <- cbind(c,d)
      }
  newdf[[i]] <- do.call(rbind, dist4)
}}


dist <- do.call(rbind, newdf)


nv <- ggplot(data=dist[dist$canwr=="nvca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("nv")
sc <- ggplot(data=dist[dist$canwr=="scnwr",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("sc")
fg <- ggplot(data=dist[dist$canwr=="fgca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("fg")
sl <- ggplot(data=dist[dist$canwr=="slnwr",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("sl")
ts <- ggplot(data=dist[dist$canwr=="tsca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("ts")
bk <- ggplot(data=dist[dist$canwr=="bkca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("bk")
cc <- ggplot(data=dist[dist$canwr=="ccnwr",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("cc")
dc <- ggplot(data=dist[dist$canwr=="dcca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("dc")
os <- ggplot(data=dist[dist$canwr=="osca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("os")
tm <- ggplot(data=dist[dist$canwr=="tmpca",])+
  geom_point(aes(x=utm_w, y=utm_n, colour=dist))+
  ggtitle("tm")

grid.arrange(nv, sc, fg, sl, ts, bk, cc, dc, os, tm, ncol=2)


hist(dist[dist$dist<=1500,]$dist, breaks=100)
