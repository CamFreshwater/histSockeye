## Study Area Figure
# March 29, 2021


library(raster)
library(dplyr)
library(sf)
library(ggplot2)


# ocean shapefile
coast_raw <- st_read(here::here("data", "shapefiles", "ne_10m_ocean.shp"))
# primary rivers shapefile
rivers_raw <- st_read(here::here("data", "shapefiles", 
                                 "ne_10m_rivers_lake_centerlines.shp"))
# secondary rivers shapefile
rivers2_raw <- st_read(here::here("data", "shapefiles", 
                                  "ne_10m_rivers_north_america.shp"))

# bounding box
min_lon <- -133.25
min_lat <- 51.5
max_lon <- -123
max_lat <- 57

coast <- coast_raw %>%
  #apply buffer to fix small errors in dataset
  st_buffer(., 0) %>%
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat)

rivers1 <- rivers_raw  %>%
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
  #convert to UTM to add buffer
  st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))

rivers2 <- rivers2_raw  %>%
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
  #convert to UTM to add buffer
  st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))


# combine primary and secondary rivers, dropping unshared columns (2 has extras)
rivers <- rbind(rivers1, rivers2[names(rivers1)]) 
rivers_sub <- rivers %>% filter(name == "Nass")

# set attribute-geometry relationship to constant to avoid errors when cropping
st_agr(rivers) = "constant"


combined_plotting <- st_union(
  coast, 
  st_buffer(rivers, 400) %>% 
    st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))
)
combined_plotting_sub <- st_union(
  coast, 
  st_buffer(rivers_sub, 200) %>% 
    st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))
)

ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_sf(data = rivers_sub, color = "black", fill = "black", size = 1) +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "darkgrey")) +
  #removes extra border
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 
