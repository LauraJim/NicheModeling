
# Function "sam.polyM" takes a random sample of environmental combinations inside
# the study area

# FUNCTIONS -------------------------------------------------------------

# Get a random sample of points inside the polygon that delimits M and extract their environmental values
## M.shp -- should be a shapefile containing the polygon that defines the area of study
## N -- sample size
## bios -- must be a stack of environmental raster layers
sam.polyM <- function(M.shp,N,bios){
  # crop and mask the environmental layers with the M polygon
  clip.M <- mask(crop(bios,M.shp),M.shp)
  # get ride of cells with NA values
  ind <- which(!is.na(clip.M[[1]][]))
  # get a random sample of indices
  sam <- sample(ind,N,replace = T)
  # choose the points corresponding to the selected indices
  Mpnts <- clip.M[][sam,]
  return(Mpnts)
}

# MAIN --------------------------------------------------------------------------------
# Laura Jimenez / # Last modification: November 2020

# Calling packages
library(raster)
library(rgdal)
library(rgeos)
library(tools)


# Read datasets ---------------------------------------------------------------

# Read environmental layers, select two
# Annual mean temperature (?C x 10)
bio1 <- raster("./ClimateData10min/bio1WH.asc")
# Annual Precipitation (mm)
bio12 <- raster("./ClimateData10min/bio12WH.asc")
# Create a single raster with as meny layers as environmental variables
stck_bios <- stack(bio1, bio12) ###

# Read M polygon
M.shp <- readOGR("./Shapefiles","C_nimbice") ###
# Get a random sample of points in M and extract its corresponding environmental values
sam.Mpnts <- sam.polyM(M.shp = M.shp,N = N,bios = stck_bios)

# Plot sampled points with occurrence points on top
