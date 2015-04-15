
# =========================================================================================
# This function takes as arguments x- and y-coordinates as vectors of length 3 which define
# the coordinates of a triangle, and returns the 3 interior angles of the triangle.
# =========================================================================================
angles.tri <- function(x,y){
  pi <- 3.14159
  ord <- order(x)
  y <- y[ord]; x <- sort(x)
  ratdiff1 <- (y[2]-y[1])/(x[2]-x[1])
  ratdiff2 <- (y[3]-y[1])/(x[3]-x[1])
  ratdiff3 <- (y[3]-y[2])/(x[3]-x[2])
  theta <- rep(0,3)
  theta[1] <- abs(atan(ratdiff1)-atan(ratdiff2))*180/pi
  theta[3] <- abs(atan(ratdiff2)-atan(ratdiff3))*180/pi
  theta[2] <- 180-abs(atan(ratdiff1)-atan(ratdiff3))*180/pi
  theta[order(ord)]
}

# =========================================================================================
# This function takes as arguments x- and y-coordinates as vectors of length 3 which define
# the coordinates of a triangle, and returns the unique point which is equidistant to all
# three points of the triangle.
# =========================================================================================

mid.tri <- function(x,y){
  y <- y[order(x)]; x <- sort(x)
  m1 <- (x[1]-x[2])/(y[2]-y[1])
  m2 <- (x[1]-x[3])/(y[3]-y[1])
  b1 <- ((y[2]^2-y[1]^2)-(x[1]^2-x[2]^2))/(2*(y[2]-y[1]))
  b2 <- ((y[3]^2-y[1]^2)-(x[1]^2-x[3]^2))/(2*(y[3]-y[1]))
  if (y[1]==y[2]) mid <- c((x[1]+x[2])/2,m2*(x[1]+x[2])/2+b2)
  if (y[1]==y[3]) mid <- c((x[1]+x[3])/2,m1*(x[1]+x[3])/2+b1)
  if ((y[1]!=y[2]) && (y[1]!=y[3]))
    mid <- solve(matrix(c(-m1,1,-m2,1),nrow=2,byrow=T))%*%matrix(c(b1,b2),nrow=2)
  mid
}

# =========================================================================================
# This function takes as arguments x- and y-coordinates as vectors of length 3 which define
# the coordinates of a triangle on the edge of the region, and returns the two intersection
# nodes from the equidistant point outside the region to the interior.
# =========================================================================================
edgenodes.tri <- function(x,y,mid,an){
  pi <- 3.14159
  obpt <- order(an)[3]
  acpts <- sort(order(an)[1:2])
  xmids <- (x[acpts]+x[obpt])/2
  ymids <- (y[acpts]+y[obpt])/2
  m1 <- (mid[2]-ymids)/(mid[1]-xmids)
  m2 <- (y[acpts[2]]-y[acpts[1]])/(x[acpts[2]]-x[acpts[1]])
  b1 <- mid[2] - m1*mid[1]
  b2 <- y[acpts[1]] - m2*x[acpts[1]]
  if (x[acpts[1]]==x[acpts[2]]){
     node1 <- c(x[acpts[1]],m1[1]*x[acpts[1]]+b1[1])
     node2 <- c(x[acpts[1]],m1[2]*x[acpts[1]]+b1[2])
  }
  if (xmids[1]==mid[1]) node1 <- c(mid[1],m2*mid[1]+b2)
  if (abs(xmids[2]-mid[1])<10^(-10)) node2 <- c(mid[1],m2*mid[1]+b2)
  if ((x[acpts[1]]!=x[acpts[2]]) && (xmids[1]!=mid[1])) node1 <-
    solve(matrix(c(-m1[1],1,-m2,1),nrow=2,byrow=T))%*%matrix(c(b1[1],b2),nrow=2)
  if ((x[acpts[1]]!=x[acpts[2]]) && (abs(xmids[2]-mid[1])>=10^(-10))) node2 <-
    solve(matrix(c(-m1[2],1,-m2,1),nrow=2,byrow=T))%*%matrix(c(b1[2],b2),nrow=2)
  matrix(c(node1,node2),nrow=2)
}
