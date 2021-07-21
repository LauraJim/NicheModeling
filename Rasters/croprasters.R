library(rgdal)
library(raster)
# read rasters with suitability index (from Niche.G), a polygon of the study area and crop
cn.wn <- raster("./Results/Catasticta_nimbice_wn_map.tif")
cn.maha <- raster("./Results/Catasticta_nimbice_maha_map.tif")
thr.wn <- raster("./Results/Threnetes_ruckeri_wn_map.tif")
thr.maha <- raster("./Results/Threnetes_ruckeri_maha_map.tif")

# read raster with environmental data
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio6 <- raster("./ClimateData10min/Bio6_america.tif")
bio12 <- raster("./ClimateData10min/bio12WH.asc")


## cropping with polygon or extents

# read polygon of study area
cn.shp <- readOGR("./Shapefiles","nimbice3")
thr.shp <- readOGR("./Shapefiles","Threnetes_ruckeri")

# crop the area to fit the polygon boundaries
cn.wnc <- mask(crop(cn.wn, cn.shp),cn.shp)
cn.mahac <- mask(crop(cn.maha, cn.shp), cn.shp)
thr.wnc <- mask(crop(thr.wn, thr.shp), thr.shp)
thr.mahac <- mask(crop(thr.maha, thr.shp), thr.shp)

writeRaster(cn.wnc,"./Rasters/Catasticta_nimbice_wn_cropped.tif", overwrite = T)
writeRaster(cn.mahac,"./Rasters/Catasticta_nimbice_maha_cropped.tif", overwrite = T)
writeRaster(thr.wnc,"./Rasters/Threnetes_ruckeri_wn_cropped.tif", overwrite = T)
writeRaster(thr.mahac,"./Rasters/Threnetes_ruckeri_maha_cropped.tif", overwrite = T)


# emap <- extent(-130, -70, 0, 40)
# bio1c <- crop(bio1, emap)
# bio12c <- crop(bio12, emap)
# cn.mahac <- crop(cn.maha, emap)
# cn.maharasc <- stack(cn.mahac, bio1c, bio12c)

bio1cn <- mask(crop(bio1, cn.shp), cn.shp)
bio6cn <- mask(crop(bio6, cn.shp), cn.shp)
bio12cn <- mask(crop(bio12, cn.shp), cn.shp)

# the extent of bio6 is off, therefore it needs to be adjusted 
# via https://stackoverflow.com/a/53440900
bio6cn.ext <- raster(vals=values(bio6cn),ext=extent(bio1cn),crs=crs(bio1cn),
              nrows=dim(bio1cn)[1],ncols=dim(bio1cn)[2])

estack <- stack(bio1cn, bio6cn.ext, bio12cn)


bio1thr <- mask(crop(bio1, thr.shp), thr.shp)
bio12thr <- mask(crop(bio12, thr.shp), thr.shp)

writeRaster(bio1cn,"./Rasters/Catasticta_nimbice_bio1_cropped.tif", overwrite = T)
writeRaster(bio12cn,"./Rasters/Catasticta_nimbice_bio12_cropped.tif", overwrite = T)
writeRaster(bio6cn.ext,"./Rasters/Catasticta_nimbice_bio6_cropped.tif", overwrite = T)

writeRaster(bio1thr,"./Rasters/Threnetes_ruckeri_bio1_cropped.tif", overwrite = T)
writeRaster(bio12thr,"./Rasters/Threnetes_ruckeri_bio12_cropped.tif", overwrite = T)

# create rasterstacks of cropped environmental rasters and models
cn.wncE <- stack(cn.wnc, bio1cn, bio12cn)
cn.mahacE <- stack(cn.mahac, bio1cn, bio12cn)
thr.wncE <- stack(thr.wnc, bio1thr, bio12thr)
thr.mahacE <- stack(thr.mahac, bio1thr, bio12thr)
