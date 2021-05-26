# Plot study area and presence points in G-space -------------------------------------------
### Plot the geographical and environmental spaces with the occurrence of the species on top
library(rgdal)
library(raster)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

occs <- data.frame(Longitude=data1$Longitud,Latitude=data1$Latitud)

# Read M polygon
M.shp <- readOGR("./Nf_models/Ecorregiones_Calibracion","CalibrationArea")
M <- spTransform(M.shp, CRS("+proj=longlat +datum=WGS84"))
M <- fortify(M)
Mcol <- rgb(0.35,0,0.2,alpha = 0.6)

# Extent of map
extent(M.shp)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## Geographic Space:
x11()
ggplot(data = world) +
  geom_sf( ) +
  theme_bw() +
  geom_map(map=M, data=M, aes(map_id=id), color="royalblue3", alpha=0.1) +
  geom_point(data = occs, aes(x = Longitude, y = Latitude), size = 2, 
             shape = 23, fill = "tomato") +
  #  annotation_scale(location = "br", width_hint = 0.2) +
  #  annotation_north_arrow(location = "br", which_north = "true", 
  #                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-10,220), ylim = c(-50,50)) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "aliceblue"))

ggsave('./Nf_models/map_Aplaci_M3_occs.png',  width = 24, height = 12, units = "cm",
       dpi = 600, pointsize = 6)
