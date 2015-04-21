
# wrote my own function
distance <- function(a,b,c,d){

  sqrt(((a-c)^2)+((b-d)^2))

}

# but here is how he did it
x <- c(1,3,1,4,5)
y <- c(5,4,3,5,1)
dat <- cbind(x,y)

dist(dat)

# trying to find the slope for kriging

a <- c(1.5,2.5,3.5,4.5,5.5)
b <- c(12.5,4.167,31.25,81.25,112.5)
dat2 <- as.data.frame(cbind(a,b))

model <- lm(dat=dat2, b ~ a-1)
