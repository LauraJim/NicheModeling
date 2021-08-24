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


## read library ----------------
library(dismo)


# End