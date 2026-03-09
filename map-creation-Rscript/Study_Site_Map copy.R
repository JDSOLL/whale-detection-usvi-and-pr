###############################################
## Clean map of the USVI and PR
###############################################

##### Load libraries
library(sf)
library(tidyverse)
library(ggplot2)
library(units)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(tigris)
library(marmap)
library(ggspatial)   # for scale bar and north arrow
library(ggrepel)

options(tigris_use_cache = TRUE)

##### Load US territories (Puerto Rico & USVI)
territories <- tigris::states(cb = TRUE, year = 2022, class = "sf")
pr_usvi <- territories[territories$STUSPS %in% c("PR", "VI"), ]

## British Virgin Islands with geodata
library(geodata)
bvi_spat  <- geodata::gadm("VGB", level = 0, path = tempdir())  # British Virgin Islands
bvi  <- sf::st_as_sf(bvi_spat)  |> sf::st_transform(4326)

##### Get bathymetry data
bathy <- getNOAA.bathy(lon1 = -68, lon2 = -63.8, lat1 = 17, lat2 = 19, resolution = 1)
bathy_df <- fortify.bathy(bathy) %>% mutate(depth_pos = pmax(0, -z))

# # Convert bathymetry to data frame for ggplot
# bathy_df <- fortify.bathy(bathy)

##### Island label coordinates (where text should be placed)
island_labels <- data.frame(
  name = c("Puerto Rico", "St. Thomas", "St. John", "St. Croix"),
  lon  = c(-66.5, -65.5, -64.5, -64.75),  # label positions
  lat  = c(18.25, 18.7, 18.7, 17.45)
)

##### Anchor points (on the islands themselves)
island_points <- data.frame(
  name = c("Puerto Rico", "St. Thomas", "St. John", "St. Croix"),
  lon  = c(-67, -64.95, -64.72, -64.75),  # approximate island centers
  lat  = c(18.2, 18.35, 18.33, 17.73)
)

##### Join label positions with anchor points
label_data <- island_labels %>%
  left_join(island_points, by = "name", suffix = c("_lab", "_anc"))

##### Make the PR line blend with the island and the rest stay black for aesthetic purposes
label_data$line_color <- c("gray90", "black", "black", "black")


##### Final plot creation
p <- ggplot() +
  # Bathymetry raster
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = depth_pos)) +
  scale_fill_gradientn(
    colours = c("lightskyblue2","dodgerblue1","royalblue2","royalblue3","royalblue4"),
    limits  = c(0, 6000),
    breaks  = seq(0, 6000, 2000),
    labels  = scales::label_number(accuracy = 1, big.mark = ""),  # <- no space/comma
    name    = "Depth (m)",
    oob     = scales::squish
  ) +
  guides(fill = guide_colorbar(reverse = TRUE)) +
  # Bathymetry contours
  geom_contour(
    data = bathy_df, aes(x = x, y = y, z = z),
    breaks = c(-50, -200, -1000),
    color = "gray20", size = 0.3
  ) +
  # Land features
  geom_sf(data = pr_usvi, fill = "gray90", color = "black") +
  geom_sf(data = bvi,     fill = "gray90", color = "black") +

  # Labels with leader lines from anchor points
  geom_label_repel(
    data = label_data,
    aes(x = lon_anc, y = lat_anc, label = name),   # <- use anchors here
    fontface = "bold", size = 4, color = "black",
    box.padding = 0.3,
    nudge_x = label_data$lon_lab - label_data$lon_anc,  # push toward label coords
    nudge_y = label_data$lat_lab - label_data$lat_anc,
    segment.color = label_data$line_color,
    segment.size = 0.4,
    fill = "white"
  ) +
  
  # Scale bar
  annotation_scale(
    location = "bl", width_hint = 0.3,
    text_col = "gray90",
    line_col = "gray90",
    style = "ticks"
  ) +
  
  # North arrow
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      line_col = "gray90",
      fill = c("gray90", "black")
    ),
    height = unit(1, "cm"),
    width = unit(1, "cm")
  ) +
  
  coord_sf(xlim = c(-68, -63.8), ylim = c(17, 19), expand = FALSE) +
  
  labs(
    title = "Study area: Puerto Rico and the\n United States Virgin Islands",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_classic() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

print(p)

# Save high-resolution figure
ggsave("/Users/JoshSoll/Desktop/Glider_Lab/Glider_Coding/RCoding/GliderMapping/CompleteMaps/V2_studysitemap.png", plot = p,
       width = 12, height = 8, dpi = 600)

