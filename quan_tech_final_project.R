
dat <- read.csv("all_birds.csv")

dat4 <- dat[dat$year==2014,]

dat4 <- na.omit(dat4)

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
            c <- b[which.min(sqrt(abs(a[j,"lat"]-b[,"lat"])^2 + abs(a[j,"long"]-b[,"long"])^2)),]
            d <- a[j,]
            c$dist <- min(sqrt(abs(a[j,"lat"]-b[,"lat"])^2 + abs(a[j,"long"]-b[,"long"])^2))
            dist4[[j]] <- cbind(c,d)
      }
  newdf[[i]] <- do.call(rbind, dist4)
}}


dist <- do.call(rbind, newdf)
