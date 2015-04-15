# This function takes an nx2 matrix of (x,y) coordinates and rotates the
# coordinates by the angle theta.  The function returns the new set of
# rotated coordinates under the name assigned to the function.  So, if
# you type: "newxy <- rotateaxis(xy,30)", newxy will be a nx2 matrix of the
# rotated (x,y) coordinates (rotated by 30 degrees).

rotateaxis <- function(xy, theta)
{
        pimult <- (theta * 2 * pi)/360
        newx <- c(cos(pimult), sin(pimult))
        newy <- c( - sin(pimult), cos(pimult))
        XY <- as.matrix(xy) %*% cbind(newx, newy)
        as.data.frame(XY)
}

