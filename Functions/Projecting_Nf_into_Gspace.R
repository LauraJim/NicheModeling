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
saveM <- "./Results/Catasticta_nimbice_map"

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
saveM2 <- "./Results/Threnetes_ruckeri_map"

result2 <- niche.G(mu = center2, Sigma = boundary2, save.map = saveM2)

x11()
plot(result2)


#### Now make plot with ggplot ------------------------------------
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
occ3 <- read.csv("./Threnetes_ruckeri_occ_bios.csv",header=T)

x11()
ggplot() +
  geom_tile(data = outppd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  # coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ3,aes(x=long, y=lat), shape = 23, fill = "yellowgreen")

ggsave('./Results/Threnetes_ruckeri_ggplot.png',  width = 24, height = 12, units = "cm",
       dpi = 600, pointsize = 6)

# Test gradient color -------------------

x11()
ggplot() +
  geom_tile(data = outppd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradientn("Suitability",limits=c(0,1), na.value = NA, +
                       guide = "colourbar", aesthetics = "fill", colors) +
  # coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ3,aes(x=long, y=lat), shape = 23, fill = "yellowgreen")
