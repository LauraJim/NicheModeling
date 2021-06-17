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

# Example 1

## Use sampled points for a plot
occ <- read.csv("./Catasticta_nimbice_bios.csv",header=T) 
pal <- c("grey50", "turquoise") # defining two colors that can be called upon

# # plot
# x11()
# # display two different graphs next to each other
# par(mfrow=c(1,2)) 
# 
# # geographic part of plot -- G-Space
# #   use long and lat of random background points; pch= display as point; 
# #   col=pal calls the second color previously defined
# plot(M.shp, col=pal[1], xlab="longitude", ylab="latitude", main="Geographic Space") 
# # add points with species location; pch=19 is form of point (full circle)
# points(occ[,1], occ[,2], pch=19, col=pal[2]) 
# 
# # add legend
# legend(x= "bottomleft",
#        legend = c("Study Area", "Occurences"),
#        lty = c(1, 0),
#        pch = c(NA, 19),
#        col = c(pal[1], pal[2]),
#        bty = "n")


# environmental part of the plot -- E-Space

# Calculate and draw the kernel using the points inside M (define M carefully)
# kernel for contour plot
fhat.M <- ks::kde(x=cbind(sam.Mpnts[,1], sam.Mpnts[,2]))
### PLOT 2: Probability of selecting a point, given the fundamental niche
# plot estimated kernel of M with points inside M
lvls1 <- c(5,25,50,75,95)
library(scales)
col.M2 <- alpha("orangered2",0.1)
M.cols <- colorRampPalette(c(alpha("white",0.1),alpha("cadetblue4",0.9)))
x11()
plot(fhat.M,display="filled.contour",cont=lvls1,main="",xlab="Mean annual temperature (°C*10)",
     ylab="Annual precipitation (mm)",col=M.cols(length(lvls1)+1))
# add points used for kernel estimation
points(sam.Mpnts[,1], sam.Mpnts[,2],col=col.M2,pch=19,cex=0.6)

#   use random points with environmental data for the plot
plot(sam.Mpnts[,1], sam.Mpnts[,2], pch=".", col=pal[1], 
     xlab="Mean Annual Temperature", ylab="Accumulated Precipitation", 
     main="Environmental Space") 
# add points of species environmental data
points(occ$bio1, occ$bio12, pch=19, col=pal[2])     

# add legend
legend(x= "topleft",
       legend = c("Study Area", "Occurences"),
       pch = c(20, 19),
       col = c(pal[1], pal[2]),
       bty = "n")



# Example 2:

# Read M polygon
M.shp <- readOGR("./threnetes_shp","Threnetes_ruckeri") ###

## Get a random sample of points in M and extract its corresponding environmental values
sam.Mpnts <- sam.polyM(M.shp = M.shp,N = 10000,bios = stck_bios)


## Use sampled points for a plot
occ <- read.csv("./Threnetes_ruckeri_occ_bios.csv",header=T) 
pal <- c("grey50", "turquoise") # defining two colors that can be called upon


# plot
x11()
par(mfrow=c(1,2))

# geographic part of plot -- G-Space
#   use long and lat of random background points; pch= display as point; 
#   col=pal calls the second color previously defined
plot(M.shp, col=pal[1], xlab="longitude", ylab="latitude", main="Geographic Space") 
# add points with species location; pch=19 is form of point (full circle)
points(occ[,1], occ[,2], pch=19, col=pal[2]) 

# add legend
legend(x= "bottomleft",
       legend = c("Study Area", "Occurences"),
       lty = c(1, 0),
       pch = c(NA, 19),
       col = c(pal[1], pal[2]),
       bty = "n")


# environmental part of the plot -- E-Space
plot(sam.Mpnts[,1], sam.Mpnts[,2], pch=".", col=pal[1], 
     xlab="Mean Annual Temperature", ylab="Accumulated Precipitation", 
     main="Environmental Space") 
points(occ$bio1, occ$bio12, pch=19, col=pal[2])  

# add legend
legend(x= "topleft",
       legend = c("Study Area", "Occurences"),
       pch = c(20, 19),
       col = c(pal[1], pal[2]),
       bty = "n")