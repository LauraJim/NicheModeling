# function E.ellipse
# 2021-05-20

# Description: "E.ellipse" ---------------
# The function E.ellipse creates confidence regions, ellipses based on the 
# Mahalanobis distance that can be used as borders for suitable environments of 
# a species.
# Various alpha levels can be chosen.

## Parameters: 
# mu = the mean of the columns that contain environmental data, such as 
#       temperature and precipitation 
# Sigma = the covariance of the environmental data linked with a species' 
#         occurrence
# Enames = character vector with the names of the environmental variables
# alpha = confidence level

## Output
# A plot with the environmental space of a species and three ellipses

# The functions code: E.ellipse --------------

E.ellipse2d <- function(Eoccs, mu, Sigma, alphas = 0.95, Enames) {
  # create a list to save different ellipses
  els <- list()
  # number of confidence levles
  la <- length(alphas)
  
  for(i in 1:la){
    # calculate the ellipses and save them into the 'els' list
    # each ellipse represents a different alpha level
    els[[i]] <- ellipse::ellipse(centre = mu, x= Sigma, level = alphas[i]) 
  }
  # create a scale of grays to color the ellipses
  pal <- gray(0:(la - 1)/la)
  
  # define values for xlim and ylim to adjust the margins of the plot
  s <- 4
  xs <- c(mu[1]-s*sqrt(Sigma[1,1]), mu[1]+s*sqrt(Sigma[1,1]))
  ys <- c(mu[2]-s*sqrt(Sigma[2,2]), mu[2]+s*sqrt(Sigma[2,2]))
  
  # create a plot that first shows the occurrences in the E-space 
  plot(Eoccs[,1], Eoccs[,2], pch=20, col= "turquoise", xlab=Enames[1],
       ylab=Enames[2], main="Environmental Space", xlim = xs, ylim = ys)
  
  # then, plot the ellipses using the different gray shades as colors
  for(i in 1:la){
    lines(els[[i]], col= pal[i], lwd = 2)
  }
  # add a legend
  legend("bottomleft",legend=c("Confidence level",alphas),lwd=2,
         col=c("white",pal),bty="n")
}

# Equivalent function for a 3d environmental space
E.ellipse3d <- function(Eoccs, mu, Sigma, alpha = 0.95){
  # define objects that contain the points on the surface of each ellipse
  elli <- ellipse3d(centre = mu, x= Sigma, level = alpha) 
  
  # initialize plotting window
  options(rgl.printRglwidget = TRUE)
  rgl.open()
  rgl.spheres(x=Eoccs[,1], y=Eoccs[,2], z=Eoccs[,3], r = 0.2, color = "#D95F02") 
  # rgl_add_axes(x, y, z, show.bbox = TRUE)
  # # Compute and draw the ellipse of concentration
  # ellips <- ellipse3d(cov(cbind(x,y,z)), 
  #                     centre=c(mean(x), mean(y), mean(z)), level = 0.95)
  shade3d(elli, col = "#D95F02", alpha = 0.1, lit = FALSE)
  # wire3d(ellips, col = "#D95F02",  lit = FALSE)
  # aspect3d(1,1,1)
  # # create a plot that shows the occurences in the environmental space 
  # plot(occ[,1], occ[,2], pch=20, col= "turquoise", xlab=Enames[1],
  #      ylab=Enames[2], main="Environmental Space", xlim = xs, ylim = ys)
  # 
  # # create a loop to write ellipse-lines with different gray colors
  # for(i in 1:la){
  #   lines(els[[i]], col= pal[i], lwd = 2)
  # }
  
}

# Main: How to use "E.ellipse" --------------

## packages needed: ellipse
library(rgl)

## Example 1: three different alpha levels

# read table with occurrences
occ <- read.csv ("./Threnetes_ruckeri_occ_GE.csv",header=T)[,-(1:2)]

# choose alpha level
alpha1 <- c(0.75, 0.9, 0.95)

# calculate the parameters
# mu calculates the means of the columns that contain the occurrences
mu1 <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
Sigma1 <- cov(occ)
# Define names for the environmental type
names1 <- c("Annual mean temperature (°C x 10)","Annual Precipitation (mm)") 

# apply function
x11()
E.ellipse2d(Eoccs=occ, mu= mu1, Sigma= Sigma1, alphas= alpha1, Enames = names1)

## Example 2:

occ <- read.csv("./Catasticta_nimbice_occ_GE.csv",header=T)[,-(1:2)]
occ1<-cbind(occ,rnorm(nrow(occ),sd=2))
# alpha-level as a sequence from 0 to 1, every 0.1 steps
alpha2 <- seq(0,1,by = 0.1)

mu2 <- colMeans(occ1)
Sigma2 <- cov(occ1)
names2 <- c("Annual mean temperature (°C x 10)","Annual Precipitation (mm)") 


E.ellipse3d(Eoccs=occ1, mu = mu2, Sigma = Sigma2)


# END

