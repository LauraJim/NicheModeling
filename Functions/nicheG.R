# Laura Jimenez and Carola Franzen
# First version: June 2020
# Last version: June 2021
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

## Output:
# The function will produce a geographical map that represents areas that have 
# suitable environmental conditions for a species. Those potential ecological 
# niche regions are colored by different degrees of suitability. The map will
# automatically be saved as a tiff and a asci file.


# the function's code: niche.G --------------------------------------

niche.G <- function(Estck, mu, Sigma) {
  
  # calculate suitability for each cell
  sui.fun <- function(cell){exp(-mahalanobis(x= c(cell[1],cell[2]), 
                                             center= mu, cov= Sigma)/2)}
  # apply this function to the whole raster layer
  suit.rast <- calc(Estck,fun=sui.fun)
  
  return(suit.rast)
}


  
# Main: How to use "niche.G" --------------- 

## libraries:

library(raster)
library(sp)
library(ggplot2)
library(ggpubr)
# needs package mvtnorm to be installed

## Read environmental layers cropped to the area of interest and stack them
bio1 <- raster(".\\ClimateData10min\\bio1WH.asc")
bio2 <- raster(".\\ClimateData10min\\bio12WH.asc")
bios <- stack(bio1,bio2)


## Example 1 - calculating the Mahalanobis distance model from species'  -------
#               occurrence with environmental data 
# Catasticta nimbice

occ2 <- read.csv("./Catasticta_nimbice_occ_GE.csv",header=T)
occ <- occ2[,-(1:2)]
  
# Set the values of the MLEs (Maximum Likelihood Estimate)
center <- colMeans(occ)
# Sigma calculates the covariance of the occurrences
boundary <- cov(occ)

# apply the function
cn.result <- niche.G(Estck = bios, mu = center, Sigma = boundary)

# write output as TIF or ASCI file
writeRaster(cn.result,"./Results/Catasticta_nimbice_maha_map.tif", overwrite = T)
writeRaster(cn.result, "./Results/Catasticta_nimbice_maha_map.asc", overwrite = T)

# plot the output
x11()
plot(cn.result)

# Another way to use the function
# Use estimated values of mu and Sigma directly
Cn.wn <- niche.G(Estck = bios, mu = c(166.1297518,1265.130825), 
                 Sigma = matrix(c(1427.054608, 7687.724366, 7687.724366, 332940.0973),ncol=2))
writeRaster(cn.result,"./Results/Catasticta_nimbice_wn_map.tif", overwrite = T)
writeRaster(cn.result, "./Results/Catasticta_nimbice_wn_map.asc", overwrite = T)


## plot in ggplot

# Read raster with output from weighted model
outp <- raster("./Results/Catasticta_nimbice_wn_map.tif")
# emap <- extent(-170, 179, -60, 80) # whole world
# emap <- extent(-140, -110, 30, 65) # USA-CAN
# emap <- extent(70, 150, 10, 55) # Asia
# emap <- extent(-15, 30, 35, 60) # Europe

# outp1 <- crop(outp, emap)
outpp <- rasterToPoints(outp)
outppd <- data.frame(outpp)
colnames(outppd) <- c("Longitude","Latitude","Suitability")


x11()
ggplot() +
  geom_tile(data = outppd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  # coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ2,aes(x=occ2[,1], y=occ2[,2]), shape = 23, fill = "yellowgreen")

ggsave('./Results/Catasticta_nimbice_nicheG_ggplot.png',  width = 24, height = 24, 
       units = "cm", dpi = 600, pointsize = 6)


## Example 2 - use values of  weighted normal distribution and Mahalanobis ----
#              distance models, calculated from fitNiche function
# Threnetes ruckeri

# read csv-file with values from both models 
thr.fitN <- read.csv("./Results/Threnetes_ruckeri_Estimated_parameters.csv",header=T)[,-1]

# use values of weighted normal distribution model
Thr.wn <- niche.G(Estck = bios, mu = thr.fitN$wn.mu, 
                  Sigma = cbind(thr.fitN$wn.sigma1, thr.fitN$wn.sigma2))

# write output as TIF or ASCI file
writeRaster(Thr.wn,"./Results/Threnetes_ruckeri_wn_map.tif", overwrite = T)
writeRaster(Thr.wn, "./Results/Threnetes_ruckeri_wn_map.asc", overwrite = T)

# use values of Mahalanobis distance model
Thr.maha <- niche.G(Estck = bios, mu = thr.fitN$maha.mu, 
                    Sigma = cbind(thr.fitN$maha.sigma1, thr.fitN$maha.sigma2))

# write output as TIF or ASCI file
writeRaster(Thr.maha,"./Results/Threnetes_ruckeri_maha_map.tif", overwrite = T)
writeRaster(Thr.maha, "./Results/Threnetes_ruckeri_maha_map.asc", overwrite = T)

# plot the output
x11()
par(mfrow=c(1,2))
plot(Thr.wn)
plot(Thr.maha)


## plot both models next to each other in ggplot (Threnetes ruckeri)

# Read presence points with geographic data
species <- read.csv("./Threnetes_ruckeri_occ_GE.csv",header=T)

# plot 1, maha
maha.map <- raster("./Results/Threnetes_ruckeri_maha_map.tif")

# create points from the raster and give column names
maha.mappt <- rasterToPoints(maha.map)
outp.maha <- data.frame(maha.mappt)
colnames(outp.maha) <- c("Longitude","Latitude","Suitability")

p1 <- ggplot() +
  geom_tile(data = outp.maha,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  geom_point(data = species,aes(x=species[,1], y=species[,2]), 
             shape = 23, fill = "yellowgreen") + 
  labs (title = "Mahalanobis distance model") +
  theme(plot.title = element_text(hjust = 0.5)) 


# plot 2, wn
wn.map <- raster("./Results/Threnetes_ruckeri_wn_map.tif")


# outp1 <- crop(outp, emap)
wn.mappt <- rasterToPoints(wn.map)
outp.wn <- data.frame(wn.mappt)
colnames(outp.wn) <- c("Longitude","Latitude","Suitability")

p2 <- ggplot() +
  geom_tile(data = outp.wn,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  geom_point(data = species,aes(x=species[,1], y=species[,2]), 
             shape = 23, fill = "yellowgreen") +
  labs (title = "Weighted normal distribution model") +
  theme(plot.title = element_text(hjust = 0.5)) 

x11()
ggarrange(p1, p2, ncol = 2, nrow = 1)

ggsave('./Results/Threnetes_ruckeri_nicheG_ggplot.png',  width = 48, height = 24,
       units = "cm", dpi = 600, pointsize = 6)




# test: try to use colorspace package for gradient filling (currently problem with factor)


# End