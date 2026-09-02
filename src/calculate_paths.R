#' This file takes some locations in space, and computes the shortest distances
#' between them, in the context of complex waterways
#' AUTHOR: Cole Brookson

#' 0 [SET UP] ------------------------------------------------------------------
library(ggplot2)
future::plan(future::multisession)
source(here::here("./src/global_funs.R"))

bc_cropped <- qs2::qs_read(here::here("./data/geodata/bc-cropped.qs2"))
locs <- readr::read_csv(here::here("./data/locationdata.csv"))

study_area_map <- ggplot(bc_cropped) +
  geom_sf(fill = "#aba9a9", colour = "#aba9a9", alpha = 1)


locs <- qs2::qs_read(here::here("./data/geodata/farm-and-sampling-locs.qs2"))
locs[, c("lat", "long")] <- sf::st_coordinates(
  sf::st_transform(locs, crs = 4326)
)

# pull the locations I want with the three farms, the three sampling sites
unique(locs$name)
temp <- locs[which(locs$type == "farm"), c("name", "lat", "long")]
ggplot(data = temp, aes(x = lat, y = long))
