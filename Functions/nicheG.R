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
  sui.fun <- function(cell){
    X <- as.vector(cell, mode = "numeric")
    sui.ind <- exp(-mahalanobis(x= X, center= mu, cov= Sigma)/2)
    return(sui.ind)
  }
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
library(rgdal)
# needs package mvtnorm to be installed

## Read environmental layers cropped to the area of interest and stack them
bio1 <- raster(".\\ClimateData10min\\bio1WH.asc")
bio2 <- raster(".\\ClimateData10min\\bio12WH.asc")
bios <- stack(bio1,bio2)


## Example 1 - calculating the Mahalanobis distance model from species'  -------
#               occurrence with environmental data 
# Catasticta nimbice

# Read matrix with geographical and environmental information of a species' occurrence
occ1 <- read.csv("./Catasticta_nimbice_occ_GE.csv",header=T)
occ2 <- occ1[,-(1:2)]
  
# Set the values of the MLEs (Maximum Likelihood Estimate)
center <- colMeans(occ2)
# Sigma calculates the covariance of the occurrences
boundary <- cov(occ2)

# apply the function
cn.result <- niche.G(Estck = bios, mu = center, Sigma = boundary)

# write output as TIF or ASCI file
writeRaster(cn.result,"./Results/Catasticta_nimbice_maha_map.tif", overwrite = T)
writeRaster(cn.result, "./Results/Catasticta_nimbice_maha_map.asc", overwrite = T)

# plot the output
x11()
plot(cn.result)

# Another way to use the function
# Use estimated values of mu and Sigma directly (weighted normal distribution)
cn.wn <- niche.G(Estck = bios, mu = c(166.1297518,1265.130825), 
                 Sigma = matrix(c(1427.054608, 7687.724366, 
                                  7687.724366, 332940.0973),ncol=2))
writeRaster(cn.wn,"./Results/Catasticta_nimbice_wn_map.tif", overwrite = T)
writeRaster(cn.wn, "./Results/Catasticta_nimbice_wn_map.asc", overwrite = T)


## plot in ggplot

# crop raster to rectangle of study area either with coordinates or polygon

# emap <- extent(-170, 179, -60, 80) # whole world
# emap <- extent(-140, -110, 30, 65) # USA-CAN
# emap <- extent(70, 150, 10, 55) # Asia
# emap <- extent(-15, 30, 35, 60) # Europe
# outp1 <- crop(cn.wn, emap)
cn.shp <- readOGR("./Shapefiles","nimbice3")
area <- crop(cn.wn, cn.shp)

# calculate raster to points for ggplot
areap <- rasterToPoints(area)
areapd <- data.frame(areap)
colnames(areapd) <- c("Longitude","Latitude","Suitability")


x11()
ggplot() +
  geom_tile(data = areapd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=4) +
  # coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ1,aes(x=occ1[,1], y=occ1[,2]), shape = 21, fill = "orange1", alpha = 0.7) +
  labs (title = "Catasticta nimbice occurrences in Middle America ") +
  labs(subtitle = "and suitable niches based on the weighted normal distribution") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 

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
  labs (title = "Weighted normal distribution model") +
  theme(plot.title = element_text(hjust = 0.5)) 

x11()
ggarrange(p1, p2, ncol = 2, nrow = 1)

ggsave('./Results/Threnetes_ruckeri_nicheG_ggplot.png',  width = 48, height = 24,
       units = "cm", dpi = 600, pointsize = 6)



# cropped ggplot:
  
thr.shp <- readOGR("./Shapefiles","Threnetes_ruckeri")


# create points from the raster and give column names
area.maha <- crop(maha.map, thr.shp)
maha.mappt <- rasterToPoints(area.maha)
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

# outp1 <- crop(outp, emap)
area.wn <- crop(wn.map, thr.shp)
wn.mappt <- rasterToPoints(area.wn)
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


# test: try to use colorspace package for gradient filling (currently problem with factor)

# End