# Function "get.Ecoord"
# Authors: Laura Jiménez and Carola Franzen
# 2021-05-18

# Description:"get.Ecoord" --------------

# The function "get.Ecoord" extracts environmental data (temperature,
# precipitation etc.) from raster-files and combines it with the coordinates of 
# tracked species. 

## Parameters 
# Estck = a raster stack with more than two climatic layers
# Gcoord = matrix with two columns whose rows are geographic coordinates of sites
#          in the study area
# Enames = character vector with the names of the environmental variables

## Output
# A matrix with four or more columns (depending on the number of layers in Estck)
# and column names: long, lat, Enames


# Function's code: get.Ecoord ----------------

get.Ecoord <- function(Estck,Gcoord, Enames) {
  # check the amount of layers in a raster and make sure they are the same as 
  # the names that are used
  if (length(Estck@layers)==length(Enames)) {  
    # extract climatic data from a raster
    a <- extract(Estck,Gcoord)  
    # combine columns of extracted coordinates with column of environmental data
    # and deletes rows with data points that are marked as NA
    b <- na.omit(cbind(Gcoord,a)) 
    # rename the columns 
    colnames(b) <- c(colnames(Gcoord),Enames) 
    return(b)
  }
  else { 
    # if the amount of layers is not the same as the amount of added column names,
    # this message is printed instead
    print("Raster stack and Enames does not have the same length.") 
  }
}

# Main: How to use "get.Ecoord" --------------

## libraries
library(dismo)

## read rasters with climatic variables
# mean annual temperature data
bio1 <- raster("./ClimateData10min/bio1WH.asc") 
# total annual precipitation data
bio12 <- raster("./ClimateData10min/bio12WH.asc") 
## combine rasters with environmental data into a single RasterStack (w/ two layers)
bios <- stack(bio1, bio12)
# add names for new columns that contain environmental data
names1 <- c("bio1", "bio12") 


# Example 1: extracting climatic values for occurrence points
# read csv of species' coordinates and delete first column (species name)
# to limit matrix to two columns with long and lat
occ <- read.csv("Catasticta_nimbice_occ_G.csv",header=T)[,-1] 

## apply function
f <- get.Ecoord(Estck=bios, Gcoord=occ, Enames=names1)

## save output
write.csv(f,file=paste0("./Catasticta_nimbice","_occ_GE.csv"),row.names = F)


# Example 2: extracting climatic values for points inside the study area

occ2 <- read.csv("Threnetes_ruckeri_occ_G.csv",header=T)[,-1]

f2 <- get.Ecoord(Estck=bios, Gcoord=occ, names1) 

write.csv(f2,file=paste0("./Threnetes_ruckeri_occ","_occ_GE.csv"),row.names = F)
