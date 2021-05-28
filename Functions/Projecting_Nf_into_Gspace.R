# Laura Jimenez
# First version: June 2020
# Last version: December 2020
# Project resulting ellipses back into G-space

# make function: --------------------------------------

niche.G <- function(mu, Sigma, save.map) {
  
  # Calculate suitabilities in each cell
  max.val <- mvtnorm::dmvnorm(x=mu,mean = mu, sigma = Sigma)
  # function that calculates the log(suitability)
  sui.fun <- function(cell){log(mvtnorm::dmvnorm(x=c(cell[1],cell[2]),mean = mu, sigma = Sigma))-log(max.val)}
  # apply this function to the whole raster layer
  suit.rast <- calc(bios,fun=sui.fun)
  # take exponential to go back to the original scale
  
  suit.rast1 <- calc(suit.rast,fun = exp,
                     filename= paste0(save.map, ".asc"), # paste puts strings 
                     # together, paste0 has no space inbetween
                     overwrite=T)
  # save a TIFF
  writeRaster(suit.rast1, paste0(save.map, ".tif"), overwrite = T)
  
  return(suit.rast1)
}
  
# Main: How to use "niche.G" --------------- 

## libraries:

library(raster)
# needs package mvtnorm to be installed

## Read environmental layers cropped to the area of interest
bio1 <- raster(".\\ClimateData10min\\bio1WH.asc")
bio2 <- raster(".\\ClimateData10min\\bio12WH.asc")
bios <- stack(bio1,bio2)

## Example 1

occ <- read.csv("./Catasticta_nimbice_bios.csv",header=T)[,-(1:2)]

# calculate parameters

# Set the values of the MLEs (Maximum Likelihood Estimate)
center <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
boundary <- cov(occ)
# define name for the maps
saveM <- "./Catasticta_nimbice_map"

result1 <- niche.G(mu = center, Sigma = boundary, save.map = saveM)

x11()
plot(result1)

## Example 2

occ2 <- read.csv("./Threnetes_ruckeri_occ_bios.csv",header=T)[,-(1:2)]

# calculate parameters

# Set the values of the MLEs (Maximum Likelihood Estimate)
center2 <- colMeans(occ2)
# Sigma calculates the covariance of the occurrences
boundary2 <- cov(occ2)
# define name for the maps
saveM2 <- "./Threnetes_ruckeri_map"

niche.G(mu = center2, Sigma = boundary2, save.map = saveM2)

x11()
plot(suit.rast1)