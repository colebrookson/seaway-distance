#' This file takes some locations in space, and computes the shortest distances
#' between them, in the context of complex waterways
#' AUTHOR: Cole Brookson

#' 0 [SET UP] ------------------------------------------------------------------
library(ggplot2)
future::plan(future::multisession)
source(here::here("./src/global_funs.R"))

bc_cropped <- qs2::qs_read(here::here("./data/geodata/bc-cropped.qs2"))
locs <- readr::read_csv(here::here("./data/locations.csv"))

# turn these locations into a helpful object (e.g. an sf object for our use)
locs <- sf::st_transform(
  # sf object
  sf::st_as_sf(locs, coords = c("long", "lat"), crs = 4326),
  # projection that works well for the area
  crs = "+proj=utm +zone=9 +datum=NAD83 +unit=m"
) |>
  dplyr::mutate(
    easting_x = sf::st_coordinates(geometry)[, "X"],
    northing_y = sf::st_coordinates(geometry)[, "Y"]
  )

# map the locations that we want to find distances between
ggplot(bc_cropped) +
  geom_sf(fill = "#aba9a9", colour = "#aba9a9", alpha = 1) +
  geom_point(
    data = locs,
    aes(
      x = easting_x,
      y = northing_y,
      fill = location_type
    ),
    colour = "black",
    shape = 21
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_better()

# 2 [GENERATE NETWORK] ---------------------------------------------------------

# create bounding box expanded a bit beyond samples
bbox <- sf::st_bbox(bc_cropped)

# turn bbox into a polygon
bbox_polygon <- sf::st_as_sfc(bbox)
bbox_polygon <- sf::st_sf(geometry = bbox_polygon, crs = sf::st_crs(bc_cropped))

# I want the inverse because I want that to be the recognized space instead of
# the land being the recognized space
inverse_area <- sf::st_difference(bbox_polygon, bc_cropped)
ggplot() +
  geom_sf(data = inverse_area) +
  geom_point(
    data = locs,
    aes(
      x = easting_x,
      y = northing_y,
      fill = location_type
    ),
    colour = "black",
    shape = 21
  ) +
  theme_better()

# generate a dense grid (adjust cellsize to trade off accuracy vs speed)
grid_sample <- sf::st_sample(
  inverse_area,
  # the size must be really large to make a fine grid - 1000 is good to see if
  # this code all runs, you can tweak how fine you want it later
  size = 50000,
  type = "regular"
) |>
  sf::st_as_sf()

# join the actual points of the sampling locations to the sampled grid
stopifnot(sf::st_crs(grid_sample) == sf::st_crs(locs))
all_nodes <- dplyr::bind_rows(
  grid_sample |> dplyr::rename(geometry = x),
  locs |> dplyr::select(geometry)
)
ggplot() +
  geom_sf(data = inverse_area) +
  geom_sf(data = grid_sample, colour = "purple", alpha = 0.1) +
  geom_sf(data = locs, colour = "blue", alpha = 0.3, size = 3) +
  theme_base()

# connect the grid
grid_connected <- nngeo::st_connect(all_nodes, all_nodes, k = 9)

# make the network itself
network <- sfnetworks::as_sfnetwork(grid_connected, directed = FALSE) |>
  tidygraph::activate("edges") |>
  dplyr::mutate(weight = sfnetworks::edge_length())

ggplot() +
  geom_sf(data = inverse_area) +
  geom_sf(
    data = network |>
      sfnetworks::activate("edges") |>
      sf::st_as_sf(),
    colour = "purple",
    alpha = 0.1
  ) +
  geom_point(
    data = locs,
    aes(
      x = easting_x,
      y = northing_y,
      fill = location_type
    ),
    colour = "black",
    shape = 21
  ) +
  theme_base()
