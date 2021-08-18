install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library("ggplot2")
theme_set(theme_bw())
library("sf")                 
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf()
ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$NAME)), " countries)"))
ggplot(data = world) + 
  geom_sf(color = "black", fill = "lightgreen")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97), expand = FALSE)
library(ggspatial)
ghana <- subset(world, admin == "Ghana")
(mainland <- ggplot(data = ghana) +
    geom_sf() +
    annotation_scale(location = "br", width_hint = 0.5) +
    annotation_north_arrow(location = "br", which_north = "true",
                           pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                           style = north_arrow_fancy_orienteering) +
    annotate(geom = "text",x=250000,y=500000,label="GHANA",fontface="italic",size=6)+
    annotate(geom = "text",x=270000,y=280000,label="Kumasi",fontface="italic",size=4)+
    annotate(geom = "point",x=270000,y=300000,shape=17,size=4)+
    annotate(geom = "point",x=220000,y=310000,shape=15,color="red",size=3)+
    annotate(geom = "point",x=320000,y=255000,shape=15,color="blue",size=3)+
    annotate(geom = "text",x=220000,y=330000,label="Bosome-Freho",fontface="italic",size=2)+
    annotate(geom = "text",x=320000,y=240000,label="Konongo",fontface="italic",size=2)+
    labs(x="Longitude",y="Latitude")+
    coord_sf(crs = st_crs(25000), xlim = c(0, 500000), ylim = c(0, 
                                                                       730000)))
ggsave("Ghana_Map.tiff",width = 8,height = 6)
