## Step 5 - generate and plot the study area of Lake Fuquene and location of the sediment records

# Load necessary packages. Install them first if you haven't already.
#install.packages(c("tidyverse", "sf", "terra", "ggplot2", "ggnewscale", "rnaturalearth", "rnaturalearthdata", "patchwork", "geodata", "rgeoboundaries"))

library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(geodata) # For downloading DEM data
#library(rgeoboundaries) # For high-resolution lake data
library(egg)
library(ggmap)

rm(list = ls())

# --- Step 1: Download and prepare geospatial data ---

# Get high-resolution data for Lake Fúquene itself
# The rgeoboundaries package is great for this, but if it fails, you can use a fallback.

## Read Lake Fuquene shapefile (created in ArcGis)
fuquene_geometry <- sf::st_read("data/GIS/fuquene_lake.shp")
st_crs(fuquene_geometry) = 4326
st_crs(fuquene_geometry)

#fuquene_geometry <- st_transform(fuquene_geometry, '+proj=longlat +datum=WGS84')
plot(fuquene_geometry)

# Create the sf object from the geometry and a data frame
fuquene_sf <- st_sf(data.frame(name = "Lake Fúquene"), geometry = st_sfc(fuquene_geometry$geometry), crs = 4326)
st_crs(fuquene_sf) 
plot(fuquene_sf)

# Get country boundaries for the main map
colombia_sf <- ne_countries(scale = "large", country = "Colombia", returnclass = "sf")
colombia_sf <- st_crop(colombia_sf, xmin = fuquene_bbox[1] - 0.1, ymin = fuquene_bbox[2] - 0.1, xmax = fuquene_bbox[3] + 0.1, ymax = fuquene_bbox[4] + 0.1)

# Get countries for the inset map
south_america <- ne_countries(continent = "South America", returnclass = "sf")

# Define the location of the lake for the inset plot
lake_point <- st_centroid(fuquene_sf)

# Create the inset plot
inset_map <- ggplot() +
  geom_sf(data = south_america, fill = "lightgrey", color = "white", size = 0.2) +
  geom_sf(data = filter(south_america, name == "Colombia"), fill = "darkgrey", color = "black", size = 0.4) +
  geom_sf(data = lake_point, color = "red", shape=15, size = 2) +
  #labs(title = "South America") +
  theme_article() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text = element_blank(), axis.ticks = element_blank()
  )
inset_map

# --- Step 4: Combine the plots using patchwork ---
final_map <- main_map +
  inset_element(inset_map, left = 0.7, bottom = 0.7, right = 0.99, top = 0.99, align_to = "panel")

# Print the final map
print(final_map)

# Save the map to a file
#ggsave("lake_fuquene_map.png", plot = final_map, width = 8, height = 8, units = "in", dpi = 300)

library(leaflet)
leaflet() %>%
  addTiles(group = "Open Street Maps") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") %>%
    addPolygons(data=fuquene_geometry) %>%
  # addPolygons(data=mysites, highlightOptions = highlightOptions(fillColor= mysites$Ecoregion, bringToFront = TRUE),
  #            popup = ~glue("{Ecoregion}")) %>%
  addLayersControl(baseGroups = c("Open Street Maps", "World Imagery"))


# Get the satellite map
fuquene_bbox <- c(-73.9, 5.3, -73.6, 5.6)

ggmap::register_google(key = "AIzaSyDGigvtrfdznthXGMcgGezNFfmvI42tutk")

register_stadiamaps("5af271ae-6451-4eb6-acc4-89ad2b1af310")
get_stadiamap()

aerial_map <- get_stadiamap(
  bbox = fuquene_bbox,
  zoom = 12,
  maptype = "stamen_terrain"
)

plot(aerial_map)
# aerial_map <- get_googlemap(
#   bbox = fuquene_bbox,
#   zoom = 12,
#   maptype = "satellite"
# )



# create two points for F7 and F9 locations
cores <- data.frame(
  core = c("F9-C", "F7"),
  lon = c(-73.732,-73.716),  
  lat = c(5.458,5.472))

cores_sf <- cores %>%
  st_as_sf(
    coords = c("lon", "lat"), # Specifies the longitude and latitude columns
    crs = 4326)              # Sets the CRS to WGS 84 (the standard for GPS data)
  

map <- ggmap(aerial_map) + # Add the lake layer
  geom_sf(data = fuquene_sf,
          fill = "steelblue",        # Polygon fill color
          alpha = 0.3,            # Transparency
          color = "white",        # Polygon border color
          linewidth = 1,
          # This is CRITICAL: Prevents CRS/aesthetic conflicts
          inherit.aes = FALSE) +
  geom_point(data=cores, aes(x=lon,y=lat,shape=core), size=8) +
  theme_article() +
  xlab("") + ylab("") +
  theme(legend.text = element_text(size=20),
        legend.title = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16))
map
  
final_sat_map <- map +
  inset_element(inset_map, left = 0.7, bottom = 0.7, right = 0.99, top = 0.99, align_to = "panel")
final_sat_map

final_map_elevation <- main_map + geom_point(data=cores, aes(x=lon,y=lat,colour=core), size=5) +
  inset_element(inset_map, left = 0.7, bottom = 0.7, right = 0.99, top = 0.99, align_to = "panel")
final_map_elevation

ggsave("outputs/for_paper/study_area_def_2.png",
       plot = final_sat_map,
       width = 14,
       height=12,
       units="in",
       dpi = 300)


