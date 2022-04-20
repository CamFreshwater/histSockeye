## Study Area Figure
# March 29, 2021


library(raster)
library(dplyr)
library(sf)
library(ggplot2)


## PREP SPATIAL DATA -----------------------------------------------------------

# ocean shapefile
# coast_raw <- st_read(here::here("data", "spatial", "shapefiles", "ne_10m_ocean.shp"))
# # primary rivers shapefile
# rivers_raw <- st_read(here::here("data", "spatial", "shapefiles", 
#                                  "ne_10m_rivers_lake_centerlines.shp"))
# # secondary rivers shapefile
# rivers2_raw <- st_read(here::here("data", "spatial", "shapefiles", 
#                                   "ne_10m_rivers_north_america.shp"))

# bounding box
# min_lon <- -133.25
# min_lat <- 51.5
# max_lon <- -123
# max_lat <- 57

# coast <- coast_raw %>%
#   #apply buffer to fix small errors in dataset
#   st_buffer(., 0) %>%
#   st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat)
# 
# rivers1 <- rivers_raw  %>%
#   st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
#   #convert to UTM to add buffer
#   st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))
# 
# rivers2 <- rivers2_raw  %>%
#   st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
#   #convert to UTM to add buffer
#   st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))
# 
# # save cropped versions
# saveRDS(coast, here::here("data", "spatial", "trim_coast_sf.rds"))
# saveRDS(rivers1, here::here("data", "spatial", "trim_rivers1_sf.rds"))
# saveRDS(rivers2, here::here("data", "spatial", "trim_rivers2_sf.rds"))


coast <- readRDS(here::here("data", "spatial", "trim_coast_sf.rds"))
rivers1 <- readRDS(here::here("data", "spatial", "trim_rivers1_sf.rds"))
rivers2 <- readRDS(here::here("data", "spatial", "trim_rivers2_sf.rds"))


# combine primary and secondary rivers, dropping unshared columns (2 has extras)
rivers <- rbind(rivers1, rivers2[names(rivers1)]) 
rivers_sub <- rivers %>% filter(name %in% c("Nass", "Bell-Irving"))

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


## PLOT ------------------------------------------------------------------------

min_lon <- -133.25
min_lat <- 53
max_lon <- -126
max_lat <- 57


# main map
main <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_sf(data = rivers_sub, color = "black", fill = "black", size = 1) +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "darkgrey")) +
  #removes extra border
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_sf(ylim = c(min_lat, max_lat), xlim = c(min_lon, max_lon))

# inset
w_can <- map_data("world", region = c("usa", "canada")) %>%
  fortify(.)
inset_map <- ggplot() +
  geom_polygon(data = w_can, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "darkgrey") + 
  labs(x = "", y = "") +
  geom_rect(aes(xmin = min_lon, xmax = max_lon, ymin = min_lat, ymax = max_lat),
            fill = NA, lty = 2, col = "red") +
  theme_linedraw() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "top",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  coord_map("azequalarea", orientation = c(60, -140, 20),
            ylim = c(40, 70), xlim = c(-155, -120))


png(here::here("outputs", "figs", "map.png"), 
    height = 7, width = 6, units = "in", res = 300)
cowplot::ggdraw(main) +
  cowplot::draw_plot(inset_map, x = 0.75, y = 0.095, vjust = -0.2, hjust = 0.05,
                     width = 0.25, height = 0.25)
dev.off()
