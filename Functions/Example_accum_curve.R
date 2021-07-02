# Laura Jim?nez
# April, 2021

# Working directory and libraries ----------------------------------------------
library(dismo)
library(maptools)
# NOTE: package 'maptools' is also needed (used inside function)

# Read functions that create accumulation curves of occurrences
source(".\\Functions\\Accumulation_curve_test.R")

## cropped suitability rasters
library(rgdal)
library(raster)
# read rasters with suitability index (from Niche.G), a polygon of the study area and crop
cn.wn <- raster("./Results/Catasticta_nimbice_wn_map.tif")
cn.maha <- raster("./Results/Catasticta_nimbice_maha_map.tif")
thr.wn <- raster("./Results/Threnetes_ruckeri_wn_map.tif")
thr.maha <- raster("./Results/Threnetes_ruckeri_maha_map.tif")

# read polygon of study area
cn.shp <- readOGR("./Shapefiles","nimbice3")
thr.shp <- readOGR("./Shapefiles","Threnetes_ruckeri")

# crop the area to fit the polygon boundaries
area.cnwn <- mask(crop(cn.wn, cn.shp),cn.shp)
area.cnmaha <- mask(crop(cn.maha, cn.shp), cn.shp)
area.thrwn <- mask(crop(thr.wn, thr.shp), thr.shp)
area.thrmaha <- mask(crop(thr.maha, thr.shp), thr.shp)

# Species name
spname <- "Catasticta_nimbice"

# Maxent model -----------------------------------------------------------------
# Occurrence data 
cn.occ <- read.csv(".\\Catasticta_nimbice_occ_G.csv",header=T)

# Apply evaluation method
# Mahalanobis
cnmaha.test <- accum.occ1(spname, output.mod=area.cnmaha, occ.pnts=cn.occ,
                         null.mod="hypergeom", conlev=0.95)

# Weighted normal model 
cnwn.test <- accum.occ1(spname, output.mod=area.cnwn, occ.pnts=cn.occ,
                          null.mod="hypergeom", conlev=0.95)

# Comparing the models ---------------------------------------------------------
# create a list with all the matrices using only the first two columns 
models <- list(cnmaha.test[,1:2], cnwn.test[,1:2])

# Apply the function that draws all the models in a single plot
comp.accplot(mods=models, nocc=155, ncells=7251,
             sp.name=spname, mods.names=c("Mahalanobis","Weigthed-Normal"),alpha=0.95)


#
#
# This function is needed to create a table that is an occurrence table with an 
# additional column
#
# Create matrices "output.mod" for MaxEnt and BioClim --------------------------
get.table <- function(occ,columns,cnames,fname=NULL){
  # Convert to points so we have long,lat values in the first two columns
  mat1 <- rasterToPoints(columns)
  colnames(mat1) <- cnames
  head(mat1)
  
  # Save the resulting matrix
  write.csv(mat1,file=paste0(".\\Rcode-hyperTest\\Example\\",fname,"_M.csv"),row.names = F)
  
  # Now repeat the previous steps with the occurrence points
  occ1 <- extract(columns,occ[,2:3])
  occ2 <- cbind(occ[,2:3],occ1)
  colnames(occ2) <- c("long","lat",sdm.name,"bio1","bio12")
  head(occ2)
  
  # Save the resulting matrix
  write.csv(occ2,file=paste0(".\\Rcode-hyperTest\\Example\\",fname,"_occ.csv"),row.names = F)

  print("Done!")
}



# END