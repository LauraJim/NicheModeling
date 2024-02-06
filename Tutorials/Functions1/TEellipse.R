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
E.ellipse3d <- function(Eoccs, mu, Sigma, alpha = 0.95, Enames){
  # define objects that contain the points on the surface of each ellipse
  elli <- ellipse3d(centre = mu, x= Sigma, level = alpha) 
  
  plot3d(x=Eoccs[,1], y=Eoccs[,2], z=Eoccs[,3], box = FALSE,
         xlab=Enames[1], ylab=Enames[2], zlab=Enames[3],
         type ="s", col = "darkorange4", size=1) 
  spheres3d(x=mu[1], y=mu[2], z=mu[3], radius=25)
  plot3d(elli, col = "darkorange4", alpha = 0.5, add = TRUE, type = "wire")
}



## read library ------------
library(rgl)


# END

