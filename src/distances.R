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

ggplot() +
    geom_sf(data = inverse_nootka) +
    geom_sf(data = network |>
        sfnetworks::activate("edges") |>
        sf::st_as_sf(), colour = "purple", alpha = 0.1) +
    geom_sf(data = samples_utm, colour = "blue", alpha = 0.3, size = 3) +
    geom_sf(
        data = farms, aes(), fill = "red", colour = "black",
        shape = 21, size = 3
    ) +
    theme_base()

# 3 [GET PATH LENGTHS] ---------------------------------------------------------
# you have to "snap" the real points to the network points
network_nodes <- sf::st_as_sf(network, "nodes")
samples_utm$nearest_node <- sf::st_nearest_feature(samples_utm, network_nodes)

sample_node_ids <- samples_utm$nearest_node

# compute all pairwise distances using manual helper function
progressr::with_progress({
    dist_table <- get_pairwise_network_distances(sample_node_ids, network)
})

# 4 [PLOT THE PATH LENGHTS] ----------------------------------------------------
# we want the geometries from out dist_table
# dist_table <- dist_table |>
#     dplyr::mutate(
#         path_geom = purrr::map(edge_paths, function(edges) {
#             net_slice <- network |>
#                 tidygraph::activate("edges") |>
#                 dplyr::slice(unlist(edge_paths)) |>
#                 sf::st_as_sf()

#             sf::st_combine(net_slice$x) |>
#                 sf::st_line_merge()
#         }, .progress = TRUE)
#     )
# dist_table$path_geom[[1]] == dist_table$path_geom[[2]]
# make an sf object from the path geometries for easier plotting
# clean_geoms <- purrr::map(dist_table$path_geom, 1)
# paths_sf <- sf::st_sf(
#     from = dist_table$from,
#     to = dist_table$to,
#     geometry = sf::st_sfc(clean_geoms),
#     crs = sf::st_crs(network)
# )

#' WHEN I COME BACK
#' just do it how i do it for kx where you activate the bit of the network you
#' want and then plot that, it doesn't need to be this complicated it hink

# plot all paths
ggplot() +
    # background
    geom_sf(data = nootka) +
    # all paths
    geom_sf(data = dist_table$edge_paths, color = "gray60", size = 0.3, alpha = 0.1) +
    # sample points
    geom_sf(data = samples_utm, color = "blue", size = 2) +
    # farm points
    geom_sf(
        data = farms_utm, shape = 21, fill = "red",
        color = "black", size = 3
    ) +
    theme_base()

plot_paths_from_node(
    from_node = samples_utm$nearest_node[1],
    to_nodes = samples_utm$nearest_node[5],
    background_sf = nootka,
    paths_sf = paths_sf,
    samples_sf = samples_utm,
    farms_sf = farms_utm
)
