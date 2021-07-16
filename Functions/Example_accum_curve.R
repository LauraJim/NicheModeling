# Laura Jim?nez
# April, 2021

# Working directory and libraries ----------------------------------------------
library(dismo)
library(maptools)
library(rgdal)
library(raster)
# NOTE: package 'maptools' is also needed (used inside function)

# Read functions that create accumulation curves of occurrences
source(".\\Functions\\Accumulation_curve_test.R")


## Example version 2 -----------------
# Parameters
# Occurrence data 
cn.occ <- read.csv(".\\Catasticta_nimbice_occ_G.csv",header=T)
thr.occ <- read.csv(".\\Threnetes_ruckeri_occ_G.csv",header=T)

# read rasters with suitability index (from Niche.G), a polygon of the study area and crop
cn.wnc <- raster("./Rasters/Catasticta_nimbice_wn_cropped.tif")
cn.mahac <- raster("./Rasters/Catasticta_nimbice_maha_cropped.tif")
thr.wnc <- raster("./Rasters/Threnetes_ruckeri_wn_cropped.tif")
thr.mahac <- raster("./Rasters/Threnetes_ruckeri_maha_cropped.tif")


# Species name
spname <- "Catasticta_nimbice"
spname2 <- "Threnetes_ruckeri"

# Apply evaluation method
# cn
cnmaha.test <- accum.occ2(spname, output.mod=cn.mahac, occ.pnts=cn.occ,
                          null.mod="hypergeom", conlev=0.95)
cnwn.test <- accum.occ2(spname, output.mod=cn.wnc, occ.pnts=cn.occ,
                        null.mod="hypergeom", conlev=0.95)

# thr 
thrmaha.test <- accum.occ2(spname, output.mod=thr.mahac, occ.pnts=thr.occ,
                           null.mod="hypergeom", conlev=0.95)
thrwn.test <- accum.occ2(spname, output.mod=thr.wnc, occ.pnts=thr.occ,
                        null.mod="hypergeom", conlev=0.95)


## Example version 3 --------------------
# Parameters
# same occurrence data and species name, and cropped raster as above

# create suitability stack for Catasticta nimbice
bio1cn <- raster("./Rasters/Catasticta_nimbice_bio1_cropped.tif")
bio12cn <- raster("./Rasters/Catasticta_nimbice_bio12_cropped.tif")

cn.wncE <- stack(cn.wnc, bio1cn, bio12cn)
cn.mahacE <- stack(cn.mahac, bio1cn, bio12cn)

# apply function
cnwn.test3 <- accum.occ3(sp.name = spname,G.occ = cn.occ,
                           suit.Estck = cn.wncE,null.mod="hypergeom",clev=0.95)
cnmaha.test3 <- accum.occ3(sp.name = spname,G.occ = cn.occ,
                           suit.Estck = cn.mahacE,null.mod="hypergeom",clev=0.95)

# create suitability stack for Threnetes ruckeri
bio1thr <- raster("./Rasters/Threnetes_ruckeri_bio1_cropped.tif")
bio12thr <- raster("./Rasters/Threnetes_ruckeri_bio12_cropped.tif")

thr.wncE <- stack(thr.wnc, bio1thr, bio12thr)
thr.mahacE <- stack(thr.mahac, bio1thr, bio12thr)

# apply function
thrwn.test3 <- accum.occ3(sp.name = spname2,G.occ = thr.occ,
                         suit.Estck = thr.wncE,null.mod="hypergeom",clev=0.95)
thrmaha.test3 <- accum.occ3(sp.name = spname2,G.occ = thr.occ,
                           suit.Estck = thr.mahacE,null.mod="hypergeom",clev=0.95)




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