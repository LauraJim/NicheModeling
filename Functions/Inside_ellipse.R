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

# Example ------------------

mu <- c(0.86,1)
sigma <- matrix(c(0.021,-0.008,-0.008,0.061),ncol=2,byrow=T)
# Define the matrix of points
cloud <- cbind(runif(5000,0,2),runif(5000,0,2))
# Use the function to determine which points from the matrix are inside the ellipse
check <- in.el(cloud,mu,sigma,0.95)

# Plot the points, ellipse and use different colors for the points inside/
# outside the ellipse
x11()
plot(check[,2],check[,3],pch=".",xlab="",ylab="",main="Points inside an ellipse")     
el <- ellipse::ellipse(x=sigma,centre = mu,level=0.95)
lines(el,col=2,lwd=2)
points(subset(check,check[,1]==1)[,2:3],pch=19,col=2)


# Example 2 ----------------
occ <- read.csv ("./Threnetes_ruckeri_occ_GE.csv",header=T)[,-(1:2)]

# mu calculates the means of the columns that contain the occurrences
mu <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
Sigma <- cov(occ)
# Define the matrix of points
cloud <- read.csv("./Threnetes_ruckeri_M_GE.csv",header=T)[,-(1:2)]

check <- in.el(cloud,mu,Sigma,0.95)

# plot
x11()
plot(check[,2],check[,3],pch=20,xlab="",ylab="",main="Points inside an ellipse")     
el <- ellipse::ellipse(x=Sigma,centre = mu,level=0.95)
lines(el,col=2,lwd=2)
points(subset(check,check[,1]==1)[,2:3],pch=19,col=2)


# example 3 ----------------

occ <- read.csv ("./Catasticta_nimbice_occ_GE.csv",header=T)[,-(1:2)]

# mu calculates the means of the columns that contain the occurrences
mu <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
Sigma <- cov(occ)
# Define the matrix of points
cloud <- read.csv("./Catasticta_nimbice_M_GE.csv",header=T)[,-(1:2)]

check <- in.el(cloud,mu,Sigma,0.95)

# plot
x11()
plot(check[,2],check[,3],pch=".",xlab="",ylab="",main="Points inside an ellipse")     
el <- ellipse::ellipse(x=Sigma,centre = mu,level=0.95)
lines(el,col=2,lwd=2)
points(subset(check,check[,1]==1)[,2:3],pch=19,col=2)

# 3D example --------------
library(rgl)

occ_GE3 <- read.csv ("./Catasticta_nimbice_occ_GE3.csv",header=T)[,-(1:2)]

# mu calculates the means of the columns that contain the occurrences
mu3 <- colMeans(occ_GE3)
# Sigma calculates the covariance of the occurrences
Sigma3 <- cov(occ_GE3)

bio1 <- raster("./ClimateData10min/bio1WH.asc") 
bio6 <- raster("./ClimateData10min/bio6.tif")
bio12 <- raster("./ClimateData10min/bio12WH.asc") 
## combine rasters with environmental data into a single RasterStack (w/ two layers)
bios <- stack(bio1, bio6, bio12)
# add names for new columns that contain environmental data
names1 <- c("bio1", "bio6", "bio12") 

M_G <- read.csv ("./Catasticta_nimbice_M_G.csv",header=T)

M_GE3 <- get.Ecoord(Estck=bios, Gcoord= M_G, Enames=names1)

write.csv(M_GE3,file=paste0("./Catasticta_nimbice","_M_GE3.csv"),row.names = F)

# Define the matrix of points
cloud3 <- read.csv("./Catasticta_nimbice_M_GE3.csv",header=T)[,-(1:2)]

open3d()
in.el(cloud = cloud3, centroid = mu3, sigma = Sigma3, alpha = 0.95)

# open3d()