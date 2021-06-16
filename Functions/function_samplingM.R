# Function: sam.polyM
# Laura Jimenez
# Last version: June 2021
# Sampling points for a geographical area

# Description: -----------------------------------------------------------
# The function "sam.polyM" takes a random sample of environmental combinations 
# inside the study area. The environmental data from a rasterstack is clipped by
# a shapefile that delimits the area before the sample is taken.

## Parameters:
# M.shp = a shapefile of the study area (polygon)
# N = the sample size
# bios = a rasterstack that contains at least two layers with environmental data

## Output:
# A matrix with two or more columns of samples that contain environmental data

# Function's code: -------------------------------------------------------------

# Get a random sample of points inside the polygon that delimits M (= )
# and extract their environmental values
sam.polyM <- function(M.shp,N,bios){
  # crop and mask the environmental layers with the M polygon
  clip.M <- mask(crop(bios,M.shp),M.shp)
  # get rid of cells with NA values = indices
  ind <- which(!is.na(clip.M[[1]][]))
  # get a random sample of indices
  sam <- sample(ind,N,replace = T)
  # choose the points corresponding to the selected indices
  Mpnts <- clip.M[][sam,]
  return(Mpnts)
}

# MAIN: How to use sam.polyM ---------------------------------------------------

## libraries
library(raster)
library(rgdal)
library(rgeos)
library(tools)


## Read datasets and prepare parameters

# Read environmental layers
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")

# Create a single raster with as many layers as environmental variables
stck_bios <- stack(bio1, bio12)

# Read M polygon
M.shp <- readOGR("./Shapefiles","C_nimbice")

# Get a random sample of points in M and extract its corresponding environmental values
# sam.Mpnts <- sam.polyM(M.shp = M.shp,N = N,bios = stck_bios)
sam.Mpnts <- sam.polyM(M.shp = M.shp,N = 10000,bios = stck_bios)

# Plot sampled points with occurrence points on top
# plot
occ <- read.csv("./Catasticta_nimbice_bios.csv",header=T) 

pal <- c("grey50", "turquoise") # defining two colors that can be called upon


## plot
x11()
par(mfrow=c(1,2)) # display two different graphs next to each other

# geographic part of plot -- G-Space
plot(M.shp, col=pal[1], xlab="longitude", ylab="latitude", main="Geographic Space") # use long and lat of random background points; pch= display as point; col=pal calls the second color previously defined
points(occ[,1], occ[,2], pch=19, col=pal[2]) # add points with species location; pch=19 is form of point (full circle)

# add legend
legend(x= "bottomleft",
       legend = c("Study Area", "Occurences"),
       lty = c(1, 0),
       pch = c(NA, 19),
       col = c(pal[1], pal[2]),
       bty = "n")


# environmental part of the plot -- E-Space
plot(sam.Mpnts[,1], sam.Mpnts[,2], pch=".", col=pal[1], xlab="Mean Annual Temperature", ylab="Accumulated Precipitation", main="Environmental Space") # use random points with environmental data
points(occ$bio1, occ$bio12, pch=19, col=pal[2]) # add points of species environmental data    

# add legend
legend(x= "topleft",
       legend = c("Study Area", "Occurences"),
       pch = c(20, 19),
       col = c(pal[1], pal[2]),
       bty = "n")


# Example 2:

# Read M polygon
M.shp <- readOGR("./threnetes_shp","Threnetes_ruckeri") ###
# Get a random sample of points in M and extract its corresponding environmental values
sam.Mpnts <- sam.polyM(M.shp = M.shp,N = 10000,bios = stck_bios)


# plot
occ <- read.csv("./Threnetes_ruckeri_occ_bios.csv",header=T) 

pal <- c("grey50", "turquoise") # defining two colors that can be called upon


## plot
x11()
par(mfrow=c(1,2)) # display two different graphs next to each other

# geographic part of plot -- G-Space
plot(M.shp, col=pal[1], xlab="longitude", ylab="latitude", main="Geographic Space") # use long and lat of random background points; pch= display as point; col=pal calls the second color previously defined
points(occ[,1], occ[,2], pch=19, col=pal[2]) # add points with species location; pch=19 is form of point (full circle)

# add legend
legend(x= "bottomleft",
       legend = c("Study Area", "Occurences"),
       lty = c(1, 0),
       pch = c(NA, 19),
       col = c(pal[1], pal[2]),
       bty = "n")


# environmental part of the plot -- E-Space
plot(sam.Mpnts[,1], sam.Mpnts[,2], pch=".", col=pal[1], xlab="Mean Annual Temperature", ylab="Accumulated Precipitation", main="Environmental Space") # use random points with environmental data
points(occ$bio1, occ$bio12, pch=19, col=pal[2]) # add points of species environmental data    

# add legend
legend(x= "topleft",
       legend = c("Study Area", "Occurences"),
       pch = c(20, 19),
       col = c(pal[1], pal[2]),
       bty = "n")