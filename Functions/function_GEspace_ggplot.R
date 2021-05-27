# function "GEspace"
# 2021-05-25
# Description:"GEspace" --------------
# This function plots two graphs with each, the geographical space and the 
# environmental space of a species occurrence in ggplot. 

## Parameters: -------------------


# Function's code:  -------------------

GEspace <- function(bckgrnd, occ, add.poly = NULL, save.p = NULL) {
  
  # rgb is color; alpha is a factor of transparency
  # Mcol <- rgb(0.35,0,0.2,alpha = 0.6)
  pal <- c("springgreen3", "tomato")
  
  if(nrow(bckgrnd)<1000)
    pch.shape = 20
  else
    pch.shape = 42
  
  if(is.null(add.poly)) {   
    # no polygon but random background points
    
    ## geographic part of plot -- G-Space
    
    bckgrnd1 <- data.frame(longitude = bckgrnd[, 1], latitude = bckgrnd[, 2])
    occ1 <- data.frame(longitude = occ[, 1], latitude = occ[, 2])
    data <- cbind(rbind(bckgrnd1[,1:2],occ1[,1:2]),
                  c(rep(1,nrow(bckgrnd1)),rep(2,nrow(occ1))))
    data2 <- data.frame(Longitude = data[, 1], Latitude = data[, 2], Type = data[,3])
    
    p1 <- ggplot(data2, aes(x = Longitude, y = Latitude, color = factor(Type), shape = factor(Type))) +
      geom_point() +
      scale_shape_manual(values=c(pch.shape,19), guide = FALSE) +
      scale_color_manual(values= c("1"=pal[1], "2"= pal[2]), guide = FALSE)
    
    ## Environmental Space:
    bckgrnd3 <- data.frame(Temperature = bckgrnd[, 3], Precipitation = bckgrnd[, 4])
    occ3 <- data.frame(Temperature = occ[, 3], Precipitation = occ [, 4])
    data3 <- cbind(rbind(bckgrnd3[,1:2],occ3[,1:2]),
                   c(rep(1,nrow(bckgrnd3)),rep(2,nrow(occ3))))
    data4 <- data.frame(Temperature = data3[, 1], Precipitation = data3[, 2], Type = data3[,3])
    
    p2 <- ggplot(data4, aes(x = Temperature, y = Precipitation, color = factor(Type), shape = factor(Type))) +
      geom_point() +
      scale_shape_manual(values=c(pch.shape,19), guide = FALSE) +
      scale_color_manual(name= "Data",
                         labels= c("Background", "Presence"),
                         values= c("1"=pal[1], "2"= pal[2])) +
      theme(legend.position = c(.05, .95), # for x, value of 0 puts it  to the left 
            # side, value of 1 to the right, for y, value of 0 puts it to the bottom,
            # value of 1 puts it to the top
            legend.justification = c("left", "top")) 
    
    
    if(!is.null(save.p)) {
      ggarrange(p1, p2, ncol = 2, nrow = 1)
      ggsave(save.p,  width = 24, height = 12, units = "cm",
             dpi = 600, pointsize = 6)
      
    }
    else{
      x11()
      ggarrange(p1, p2, ncol = 2, nrow = 1)
    }
  }
  
  else { 
    ## use polygon in map
    
    # geographical space:
    occ5 <- data.frame(Longitude = occ[, 1], Latitude = occ [, 2])
    M <- spTransform(add.poly, CRS("+proj=longlat +datum=WGS84"))
    # This function turns a map into a data frame that can more easily be plotted 
    # with ggplot2.
    M <- fortify(M)
    # takes id that is a "character" and converts it to a number
    M$id = as.numeric(M$id)
    
    ## Extent of map
    # where to cut the map
    Mext <- extent(add.poly)
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    class(world)
    
    
    p3 <- ggplot(data = world) +
      geom_sf( ) +
      theme_bw() +
      # alpha dictates the transparency of the polygon, 0.1 is very transparent 0.5 is medium
      geom_map(map=M, data=M, aes(map_id=id), color= pal[1], fill = pal[1], alpha = 0.4) +
      geom_point(data = occ5, aes(x = Longitude, y = Latitude), size = 2, 
                 shape = 23, fill = pal[2]) +
      #  annotation_scale(location = "br", width_hint = 0.2) +
      #  annotation_north_arrow(location = "br", which_north = "true", 
      #                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
      #                         style = north_arrow_fancy_orienteering) +
      coord_sf(xlim = c(Mext[1] - 10 ,Mext[2] + 10), ylim = c(Mext[3] - 10,Mext[4] + 10), expand = FALSE) +
      theme(panel.grid.major = element_line(color = "white"),
            panel.background = element_rect(fill = "aliceblue"))
    
    ## Environmental Space:
    bckgrnd3 <- data.frame(Temperature = bckgrnd[, 3], Precipitation = bckgrnd[, 4])
    occ3 <- data.frame(Temperature = occ[, 3], Precipitation = occ [, 4])
    data3 <- cbind(rbind(bckgrnd3[,1:2],occ3[,1:2]),
                   c(rep(1,nrow(bckgrnd3)),rep(2,nrow(occ3))))
    data4 <- data.frame(Temperature = data3[, 1], Precipitation = data3[, 2], Type = data3[,3])
    
    p4 <- ggplot(data4, aes(x = Temperature, y = Precipitation, color = factor(Type), shape = factor(Type))) +
      geom_point() +
      scale_shape_manual(values=c(pch.shape,19), guide = FALSE) +
      scale_color_manual(name= "Data",
                         labels= c("Background", "Presence"),
                         values= c("1"=pal[1], "2"= pal[2])) +
      theme(legend.position = c(.05, .95), # for x, value of 0 puts it  to the left 
            # side, value of 1 to the right, for y, value of 0 puts it to the bottom,
            # value of 1 puts it to the top
            legend.justification = c("left", "top"))
    
    
    if(!is.null(save.p)) {
      ggarrange(p3, p4, ncol = 2, nrow = 1)
      ggsave(save.p,  width = 24, height = 12, units = "cm",
             dpi = 600, pointsize = 6)
      
    }
    else{
      x11()
      ggarrange(p3, p4, ncol = 2, nrow = 1)
    }
    
  }
  
}

# Main: How to use "GEspace" function ------------------------

library(rgdal)
library(raster)
library(ggplot2)
# library(sf) # does not need to be run but installed
library(rnaturalearth) # package rgeos needs to be installed to run it (no need to run library)
library(rnaturalearthdata)
library(ggpubr)



## Example 1:

ranpoints <- read.csv("./Catasticta_nimbice_M_bios.csv",header=T)
species <- read.csv("./Catasticta_nimbice_bios.csv",header=T)
shp <- readOGR("./Shapefiles","nimbice3")

GEspace(bckgrnd=ranpoints, occ=species, add.poly=shp, save.p = "map1.png")
GEspace(bckgrnd=ranpoints, occ=species, add.poly=shp)
GEspace(bckgrnd=ranpoints, occ=species, save.p = "map2.png")
GEspace(bckgrnd=ranpoints, occ=species)

## Example 2:
ranpoints2 <- read.csv("./Threnetes_ruckeri_M_bios.csv",header=T) # use random points as background
species2 <- read.csv("./Threnetes_ruckeri_occ_bios.csv",header=T) # occurence of species with environmental data (previously created file)
shp2 <- readOGR("./Shapefiles","Threnetes_ruckeri")

GEspace(bckgrnd = ranpoints2, occ = species2, add.poly = shp2, save.p = "map3.png")
GEspace(bckgrnd = ranpoints2, occ = species2, add.poly = shp2)
GEspace(bckgrnd = ranpoints2, occ = species2, save.p = "map4.png")
GEspace(bckgrnd = ranpoints2, occ = species2)
