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

# rasters for rasterstack
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")
cn.maharas <- stack(cn.maha, bio1, bio12)

emap <- extent(-130, -70, 0, 40)
bio1c <- crop(bio1, emap)
bio12c <- crop(bio12, emap)
cn.mahac <- crop(cn.maha, emap)
cn.maharasc <- stack(cn.mahac, bio1c, bio12c)

bio1cr <- mask(crop(bio1, cn.shp), cn.shp)
bio12cr <- mask(crop(bio12, cn.shp), cn.shp)
cn.maharc <- stack(area.cnmaha, bio1cr, bio12cr)

# Species name
spname <- "Catasticta_nimbice"

# Maxent model -----------------------------------------------------------------
# Occurrence data 
cn.occ <- read.csv(".\\Catasticta_nimbice_occ_G.csv",header=T)

# Apply evaluation method
# Mahalanobis
cnmaha.test <- accum.occ1(spname, output.mod=area.cnmaha, occ.pnts=cn.occ,
                         null.mod="hypergeom", conlev=0.95)
cnmaha.test3 <- accum.occ3(sp.name = spname,G.occ = cn.occ,
                           suit.Estck = cn.maharasc,null.mod="hypergeom",conlev=0.95)

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
get.table <- function(G.occ,suit.Estck){
  # Convert to points so we have long,lat values in the first two columns
  mat1 <- rasterToPoints(suit.Estck)
  # order of index, ask order of rows order (range of indexes)
  iord <- order(mat1[,3], decreasing = T)
  # create new matrix with new order by suitability rows (high to low)
  mat2 <- mat1[iord,]
  
  # Now repeat the previous steps with the occurrence points
  occ1 <- extract(suit.Estck,G.occ[,2:3])
  occ2 <- na.omit(cbind(G.occ[,2:3],occ1))
  # sort the values of vector
  # order of index, ask order of rows order (range of indexes)
  iord2 <- order(occ2[,3], decreasing = T) 
  occ3 <- occ2[iord2,]
  colnames(occ3) <- colnames(mat2)
  mat3 <- rbind(occ3, mat2)
  mat4 <- cbind(mat3, Type=c(rep(1, nrow(occ3)), rep(0, nrow(mat2))))

  return(mat4)
}

cn.occG <- read.csv("./Catasticta_nimbice_occ_G.csv",header=T)

bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")
cn.maharas <- stack(cn.maha, bio1, bio12)

test.cnmaha <- get.table(G.occ = cn.occG, suit.Estck = cn.maharas)
dim(test.cnmaha[test.cnmaha$Type==1,])

# END