
dat <- read.csv("all_birds.csv")

dat4 <- dat[dat$year==2014,]

jdate4 <- unique(dat4$jdate)


list4 <- list()

for(i in 1:length(jdate4)){
  list4[[i]] <- dat4[dat4$jdate==jdate4[i],]
}


dist4 <- list()

for(i in 1:length(jdate4)){
  df <- as.data.frame(list4[[i]])
  u <- unique(df$obs)
  if (length(u)>1 & length(df)>1) {
    a <- df[df$obs=="A",]
    b <- df[df$obs=="N",]
            for(j in 1:nrow(a)){
            c <- b[which.min(sqrt(abs(a[j,"lat"]-b[,"lat"])^2 + abs(a[j,"long"]-b[,"long"])^2)),]
            d <- a[j,]
            c$dist <- min(sqrt(abs(a[j,"lat"]-b[,"lat"])^2 + abs(a[j,"long"]-b[,"long"])^2))
            dist4[[j]] <- cbind(c,d)
            }
      }
}


dist <- do.call(rbind, dist4)
