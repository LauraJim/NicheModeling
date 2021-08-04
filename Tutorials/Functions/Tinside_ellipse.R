# May, 2017
# Laura Jimenez

# Function that determines which points from a given matrix are inside
# an ellipse -------

## Description: --------------
# The function `in.el` creates a matrix that determines how many points are 
# within a confidence region bordered by an ellipse. The ellipse is calculated 
# based on a species' environmental data and a confidence level. The points 
# inside the ellipse are potential data points of the niche in which a species 
# may thrive.

## Parameters ----------
## cloud == matrix that contains the coordinates of the points to evaluate
## centroid == vector with the coordinates of the centroid of the ellipse
## sigma == covariance matrix that defines de ellipse
## alpha == confidence level


in.el <- function(cloud,centroid,sigma,alpha){
  # step 1: calculate de Mahalanobis distance
  maha <- mahalanobis(x=cloud,center=centroid,cov=sigma)
  # step 2: a point is inside the confidence region (1-alpha=confidence%) if
  # the distance divided by the quantile of a Chi-square variable with k d.f. is less than 1
  chi.q <- qchisq(alpha,ncol(cloud))
  check <- (maha/chi.q) <= 1
  cloud <- cbind(rep(1,nrow(cloud))*check,cloud)
    return(as.matrix(cloud))
}


## read library -------------
library(rgl)

# End