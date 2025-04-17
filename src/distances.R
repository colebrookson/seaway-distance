# distances.R
# Author: Cole Brookson
# Date: 2025-04-17
# Description: This script contains functions and calculations related to
# computing distances betwen points in seaway distance


#' 0 [SET UP] ------------------------------------------------------------------
library(ggplot2)

# get functions I'll want
source(here::here("./src/global_funs.R"))

farms <- readr::read_csv(here::here("./data/farmlocations.csv"))
samples <- readr::read_csv(here::here("./data/dispersalcoordinates.csv")) |>
    dplyr::filter(region == "Nootka")
geo_data <- sf::read_sf(here::here("./data/geodata/kxx-shapefile/"))

#' 1 [CREATE GEOMETRIES] -------------------------------------------------------
samples <- samples |>
    sf::st_as_sf(coords = c("dd_long", "dd_lat"), crs = 4326) # assuming WGS84
farms <- farms |>
    sf::st_as_sf(coords = c("dd_long", "dd_lat"), crs = 4326) # assuming WGS84

##' 1.1 [TRANSFER TO UTM] ------------------------------------------------------
# easiest to do this in UTM, so we'll change the coordinates first
samples_utm <- sf::st_transform(
    samples,
    crs = "+proj=utm +zone=9 +datum=NAD83 +unit=m"
)
farms_utm <- sf::st_transform(
    farms,
    crs = "+proj=utm +zone=9 +datum=NAD83 +unit=m"
)

##' 1.2 [EXTRACT COORDINATES AND COMBINE] --------------------------------------
coords_sample <- sf::st_coordinates(samples_utm)
samples_utm <- cbind(
    samples_utm,
    x_utm = coords_sample[, "X"],
    y_utm = coords_sample[, "Y"]
)
coords_farm <- sf::st_coordinates(farms_utm)
farms_utm <- cbind(
    farms_utm,
    x_utm = coords_farm[, "X"],
    y_utm = coords_farm[, "Y"]
)

##' 1.2 [TRIM SHAPEFILE] ---------------------------------------------------------
# need this in an sf type
geo_data_sf_bc <- sf::st_as_sf(geo_data)
ggplot() +
    geom_sf(data = geo_data_sf_bc)
# switch to UTM
nootka_utm <- sf::st_transform(geo_data_sf_bc,
    crs = "+proj=utm +zone=9 +datum=NAD83 +unit=m"
)

# get difference between top and bottom / left and right of the sampled area
x_dist <- max(samples_utm$x_utm) - min(samples_utm$x_utm)
y_dist <- max(samples_utm$y_utm) - min(samples_utm$y_utm)
nootka <- sf::st_crop(nootka_utm,
    xmin = min(samples_utm$x_utm) - 0.09 * x_dist,
    xmax = max(samples_utm$x_utm) + 0.09 * x_dist,
    ymin = min(samples_utm$y_utm) - 0.09 * y_dist,
    ymax = max(samples_utm$y_utm) + 0.09 * y_dist
)
sf::st_bbox(nootka)
ggplot() +
    geom_sf(data = nootka) +
    geom_sf(data = samples_utm, colour = "blue", alpha = 0.3, size = 3) +
    geom_sf(
        data = farms, aes(), fill = "red", colour = "black",
        shape = 21, size = 3
    ) +
    theme_base()

# 2 [GENERATE NETWORK] ---------------------------------------------------------
# create bounding box expanded a bit beyond samples
bbox <- sf::st_bbox(nootka)
# turn bbox into a polygon
bbox_polygon <- sf::st_as_sfc(bbox)
bbox_polygon <- sf::st_sf(geometry = bbox_polygon, crs = sf::st_crs(nootka))

# I want the inverse because I want that to be the recognized space instead of
# the land being the recognized space
inverse_nootka <- sf::st_difference(bbox_polygon, nootka)

ggplot() +
    geom_sf(data = inverse_nootka) +
    geom_sf(data = samples_utm, colour = "blue", alpha = 0.3, size = 3) +
    geom_sf(
        data = farms, aes(), fill = "red", colour = "black",
        shape = 21, size = 3
    ) +
    theme_base()

# generate a dense grid (adjust cellsize to trade off accuracy vs speed)
grid_sample <- sf::st_sample(
    inverse_nootka,
    # the size is really large to make a fine grid
    size = 5000, type = "regular"
) |>
    sf::st_as_sf()

ggplot() +
    geom_sf(data = inverse_nootka) +
    geom_sf(data = grid_sample, colour = "purple", alpha = 0.1) +
    geom_sf(data = samples_utm, colour = "blue", alpha = 0.3, size = 3) +
    geom_sf(
        data = farms, aes(), fill = "red", colour = "black",
        shape = 21, size = 3
    ) +
    theme_base()

# connect the grid
grid_connected <- nngeo::st_connect(grid_sample, grid_sample, k = 9)

# make the network itself
network <- sfnetworks::as_sfnetwork(grid_connected, directed = FALSE) |>
    tidygraph::activate("edges") |>
    dplyr::mutate(weight = sfnetworks::edge_length())
