# function e.ellipse
# 2021-05-20

# Description: "e.ellipse" ---------------
# The function e.ellipse creates confidence regions, ellipses based on the 
# Mahalanobis distance that can be used as borders for suitable environments of 
# a species
# Currently it only works with three alpha levels

## Parameters: 
# mu = the mean of the columns that contain environmental data, such as 
# temperature and precipitation 
# Sigma = the covariance of the occurence
# enames = character vector with the names of the environmental variables

## Output
# A plot with the environmental space of a species and three ellipses

# The functions code: e.ellipse --------------

e.ellipse <- function(mu, Sigma, alpha = 0.95, enames) {
  # create a list
  els <- list()
  
  la <- length(alpha)
  
  for(i in 1:la){
    els[[i]] <- ellipse::ellipse(centre = mu, x= Sigma, level = alpha[i]) 
    # create a loop that inputs the ellipses into the empty list based on the 
    # chosen alpha levels (loop repeats for as many alpha levels are put in)
  }
  # create a scale of grays for the ellipses
  pal <- gray(0:(la - 1)/la)
  
  # define values for xlim and ylim to adjust the margins
  xs <- c(min(occ[,1]) -20, max(occ[,1]) +20)
  ys <- c(min(occ[,2]) -200, max(occ[,2]) +200)
  
  # create a plot that shows the occurences in the environmental space 
  plot(occ[,1], occ[,2], pch=20, col= "turquoise", xlab=enames[1],
       ylab=enames[2], main="Environmental Space", xlim = xs, ylim = ys)
  
  # create a loop to write ellipse-lines with different gray colors
  for(i in 1:la){
    lines(els[[i]], col= pal[i], lwd = 2)
  }
  
}

# Main: How to use "e.ellipse" --------------

## packages needed: ellipse


## Example 1: three different alpha levels

# read table with occurrences
occ <- read.csv ("./Threnetes_ruckeri_occ_bios.csv",header=T)[,-(1:2)]

# choose alpha level
alpha <- c(0.75, 0.9, 0.95)

# calculate the parameters
# mu calculates the means of the columns that contain the occurrences
mu <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
Sigma <- cov(occ)
# Define names for the environmental type
names1 <- c("Annual mean temperature (°C x 10)","Annual Precipitation (mm)") 

# apply function  
f <- e.ellipse(mu, Sigma, alpha, enames = names1)


## Example 2:

occ <- read.csv("./Catasticta_nimbice_bios.csv",header=T)[,-(1:2)]

# alpha-level as a sequence from 0 to 1, every 0.1 steps
alpha2 <- seq(0,1,by = 0.1)

mu2 <- colMeans(occ)
Sigma2 <- cov(occ)
names2 <- c("Annual mean temperature (°C x 10)","Annual Precipitation (mm)") 

f <- e.ellipse(mu = mu2, Sigma = Sigma2, alpha = alpha2, enames = names2)
