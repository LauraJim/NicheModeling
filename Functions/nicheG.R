# Laura Jimenez
# First version: June 2020
# Last version: December 2020
# Project resulting ellipses back into G-space

# Description: -----------------------------------
# The function niche.G projects ellipses that define suitable environments for a 
# species on a map as potential niches. The regions in the geographical space
# are colored by different degrees of suitability.

## Parameters:
# Estck = a raster stack with more than two climatic layers
# mu = the mean of the columns that contain environmental data, such as 
#       temperature and precipitation 
# Sigma = the covariance of the environmental data linked with a species' 
#         occurrence
# save.map = set the location and name for saving the map

## Output:
# The function will produce a geographical map that represents areas that have 
# suitable environmental conditions for a species. Those potential ecological 
# niche regions are colored by different degrees of suitability. The map will
# automatically be saved as a tiff and a asci file.

# the function's code: niche.G --------------------------------------

niche.G <- function(Estck, mu, Sigma, save.map) {
  
  # Calculate suitabilities in each cell
  max.val <- mvtnorm::dmvnorm(x=mu,mean = mu, sigma = Sigma)
  # function that calculates the log(suitability)
  sui.fun <- function(cell){log(mvtnorm::dmvnorm(x=c(cell[1],cell[2]),mean = mu, 
                                                 sigma = Sigma))-log(max.val)}
  # apply this function to the whole raster layer
  suit.rast <- calc(Estck,fun=sui.fun)
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
library(sp)
library(ggplot2)
# needs package mvtnorm to be installed

## Read environmental layers cropped to the area of interest
bio1 <- raster(".\\ClimateData10min\\bio1WH.asc")
bio2 <- raster(".\\ClimateData10min\\bio12WH.asc")
bios <- stack(bio1,bio2)


## Example 1

occ <- read.csv("./Catasticta_nimbice_occ_GE.csv",header=T)[,-(1:2)]

# calculate parameters

# Set the values of the MLEs (Maximum Likelihood Estimate)
center <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
boundary <- cov(occ)
# define name for the maps
saveM <- "./Results/Catasticta_nimbice_map"

result1 <- niche.G(Estck = bios, mu = center, Sigma = boundary, save.map = saveM)

x11()
plot(result1)


## Example 2

occ2 <- read.csv("./Threnetes_ruckeri_occ_GE.csv",header=T)[,-(1:2)]

# calculate parameters

# Set the values of the MLEs (Maximum Likelihood Estimate)
center2 <- colMeans(occ2)
# Sigma calculates the covariance of the occurrences
boundary2 <- cov(occ2)
# define name for the maps
saveM2 <- "./Results/Threnetes_ruckeri_map"

result2 <- niche.G(Estck, mu = center2, Sigma = boundary2, save.map = saveM2)

x11()
plot(result2)


#### Now make plot with ggplot ------------------------------------

## Example 1: Threnetes ruckeri

# Read raster with output from weighted model
outp <- raster("./Results/Threnetes_ruckeri_map.tif")
# emap <- extent(-170, 179, -60, 80) # whole world
# emap <- extent(-140, -110, 30, 65) # USA-CAN
# emap <- extent(70, 150, 10, 55) # Asia
# emap <- extent(-15, 30, 35, 60) # Europe

# outp1 <- crop(outp, emap)
outpp <- rasterToPoints(outp)
outppd <- data.frame(outpp)
colnames(outppd) <- c("Longitude","Latitude","Suitability")

# Read presence points 
occ3 <- read.csv("./Threnetes_ruckeri_occ_GE.csv",header=T)

x11()
ggplot() +
  geom_tile(data = outppd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  # coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ3,aes(x=occ3[,1], y=occ3[,2]), shape = 23, fill = "yellowgreen")

ggsave('./Results/Threnetes_ruckeri_ggplot.png',  width = 24, height = 24, units = "cm",
       dpi = 600, pointsize = 6)


## Example 2: Catasticta nimbice

# Read raster with output from weighted model
outp2 <- raster("./Results/Catasticta_nimbice_map.tif")
outpp2 <- rasterToPoints(outp2)
outppd2 <- data.frame(outpp2)
colnames(outppd2) <- c("Longitude","Latitude","Suitability")

# Read presence points 
occ4 <- read.csv("./Catasticta_nimbice_occ_GE.csv",header=T)

# plot
x11()
ggplot() +
  geom_tile(data = outppd2,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  # coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ4,aes(x=occ4[,1], y=occ4[,2]), shape = 23, fill = "yellowgreen")

ggsave('./Results/Catasticta_nimbice_ggplot.png',  width = 24, height = 24, units = "cm",
       dpi = 600, pointsize = 6)



# Test gradient color - still work in progress -------------------
library(colorspace)

# issue with the discrete color filling, data needs to be converted via factor 
# but this does not work out (does not finish processing plot)

x11()
ggplot() +
  geom_tile(data = outppd2,aes(x=Longitude, y=Latitude, color= Suitability)) +
#  theme_bw() +
#  geom_density(alpha = 0.6) +
  scale_color_discrete_sequential(palette = "BluYl") + 
  # scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
  #                     mid='slateblue1', high = 'slateblue4',na.value = NA,
  #                     midpoint = 0.5, n.breaks=4) +
  geom_point(data = occ4,aes(x=occ4[,1], y=occ4[,2]), shape = 23, fill = "yellowgreen")



#
## Use function to create maps using output from weighted model and mahalanobis model
# example Catasticta nimbice

cn.fitN <- read.csv("./Results/Catasticta_nimbice_Estimated_parameters.csv",header=T)[,-1]

# weighted
# define name for the maps
saveM.wn <- "./Results/Catasticta_nimbice_wn_map"

Cn.wn <- niche.G(Estck = bios, mu = cn.fitN$wn.mu, 
                 Sigma = cbind(cn.fitN$wn.sigma1,cn.fitN$wn.sigma2),
                 save.map = saveM.wn)

# # example maha
# # define name for the maps
# saveM.maha <- "./Results/Catasticta_nimbice_maha_map"
# 
# Cn.maha <- niche.G(Estck = bios, mu = cn.fitN$maha.mu, 
#                    Sigma = cbind(cn.fitN$maha.sigma1,cn.fitN$maha.sigma2), 
#                    save.map = saveM.maha)
# 

# example Threnetes ruckeri
# weighted
cn.fitN <- read.csv("./Results/Threnetes_ruckeri_Estimated_parameters.csv",header=T)

center2.wn <- cn.fitN[1:2,1]
center2.wn2 <- sapply(center2.wn,as.numeric)

boundary2.wn <- rbind(cn.fitN[1,2:3], cn.fitN[2,2:3])
boundary2.wn2 <- sapply(boundary2.wn,as.numeric)

saveM2.wn <- "./Results/Threnetes_ruckeri_wn_map"

Thr.wn <- niche.G(Estck = bios, mu = center2.wn2, Sigma = boundary2.wn2, save.map = saveM2.wn)

# example maha
center2.maha <- cn.fitN[1:2,4]
center2.maha2 <- sapply(center2.maha,as.numeric)

boundary2.maha <- rbind(cn.fitN[1,5:6], cn.fitN[2,5:6])
boundary2.maha2 <- sapply(boundary2.maha,as.numeric)

saveM2.maha <- "./Results/Threnetes_ruckeri_maha_map"

Thr.maha <- niche.G(Estck = bios, mu = center2.maha2, Sigma = boundary2.maha2, save.map = saveM2.maha)
