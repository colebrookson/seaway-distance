#' This file takes some locations in space, and computes the shortest distances
#' between them, in the context of complex waterways
#' AUTHOR: Cole Brookson

#' 0 [SET UP] ------------------------------------------------------------------
library(ggplot2)
library(patchwork)

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
all_nodes <- grid_sample |> dplyr::rename(geometry = x)
ggplot() +
  geom_sf(data = inverse_area) +
  geom_sf(data = grid_sample, colour = "#cd99ed", alpha = 0.1, size = 0.01) +
  geom_sf(data = locs, colour = "blue", alpha = 0.3, size = 3) +
  theme_better()

# connect the grid
# note that k = 9 means a 4 orthogonal and 4 diagonal connection cells, so each
# edge is like roughly 45 degrees, so if we make it k = 25, we have a 5 x 5
# block which gives 16 possible directions, so it's less likely too do
# "inefficient" style of distances
grid_connected <- nngeo::st_connect(all_nodes, all_nodes, k = 25)

land <- sf::st_cast(sf::st_geometry(bc_cropped), "POLYGON")
len <- as.numeric(sf::st_length(grid_connected))
grid_connected <- grid_connected[len > 0]
len <- len[len > 0]

step <- min(len)
long_idx <- which(len > step * 1.5)

message(length(grid_connected), " edges, ", length(long_idx), " long")

if (length(long_idx) > 0) {
  hits <- lengths(sf::st_intersects(grid_connected[long_idx], land)) > 0
  message(sum(hits), " crossing land")
  if (any(hits)) grid_connected <- grid_connected[-long_idx[hits]]
}

stopifnot(length(grid_connected) > 0)

# make the network itself
network <- sfnetworks::as_sfnetwork(grid_connected, directed = FALSE) |>
  sfnetworks::st_network_blend(locs) |>
  tidygraph::activate("edges") |>
  dplyr::mutate(weight = sfnetworks::edge_length())
qs2::qs_save(network, here::here("network-50k-k25.qs2"))

# check that there's no land problems
crossings <- lengths(sf::st_intersects(
  sf::st_geometry(edges_sf)[unique(unlist(path_df$edge_paths))],
  land
)) >
  0

sum(crossings)

ggplot() +
  geom_sf(data = inverse_area) +
  geom_sf(
    data = network |>
      sfnetworks::activate("edges") |>
      sf::st_as_sf(),
    colour = "#cd99ed",
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
  theme_better()

# 3 [GET PATH LENGTHS] ---------------------------------------------------------
network_nodes <- sf::st_as_sf(network, "nodes")
locs$network_nodes <- match(locs$location_name, network_nodes$location_name)
stopifnot(!anyNA(locs$network_nodes))
# confirm all six are mutually reachable before computing anything
comp <- network |>
  tidygraph::activate("nodes") |>
  dplyr::mutate(comp = tidygraph::group_components()) |>
  tidygraph::pull(comp)
stopifnot(length(unique(comp[locs$network_nodes])) == 1)

# calculate total path length for each route (full factorial)
edge_weights <- network |>
  tidygraph::activate("edges") |>
  tidygraph::pull(weight)

# set source nodes as the sampling sites
node_ids <- locs$network_nodes[which(locs$location_type == "sampling")]

# choose sequential or parallel execution
use_parallel <- FALSE
map_fun <- if (use_parallel) furrr::future_map else purrr::map

# calculate all pairwise paths
all_paths <- map_fun(
  node_ids,
  function(from_node) {
    # helper function use here
    to_nodes <- get_destination_nodes(from_node, "sampling", locs)
    paths <- sfnetworks::st_network_paths(
      x = network,
      from = from_node,
      to = to_nodes,
      weights = "weight"
    )
    tibble::tibble(
      from = from_node,
      to = to_nodes,
      edge_paths = paths$edge_paths
    )
  }
)

# flatten into a single tibble
path_df <- dplyr::bind_rows(all_paths)

# exclude self paths
path_df <- path_df |>
  dplyr::filter(from != to) |>
  dplyr::mutate(n_edges = lengths(edge_paths))

if (any(path_df$n_edges == 0)) {
  warning(sum(path_df$n_edges == 0), " unreachable pair(s) dropped")
  path_df <- dplyr::filter(path_df, n_edges > 0)
}

path_df <- path_df |>
  dplyr::mutate(
    path_length = purrr::map_dbl(edge_paths, ~ sum(edge_weights[.x]))
  )

# 4 [PLOT THE PATH LENGHTS] ----------------------------------------------------
farm_cols <- setNames(
  c("#0072B2", "#D55E00", "#009E73"),
  sort(unique(farms$location_name))
)
type_cols <- c(
  "farm" = "#E8C547",
  "sampling" = "#333333",
  "focal" = "#CC3311" # red for the panel's source site
)

edges_sf <- network |>
  sfnetworks::activate("edges") |>
  sf::st_as_sf()

path_labels <- path_df |>
  dplyr::mutate(
    from_name = locs$location_name[match(from, locs$network_nodes)],
    to_name = locs$location_name[match(to, locs$network_nodes)],
    mid_edge = purrr::map_int(edge_paths, \(e) {
      as.integer(e[ceiling(length(e) / 2)])
    }),
    label = sprintf("%.1f km", path_length / 1000)
  )

path_labels <- sf::st_sf(
  path_labels,
  geometry = sf::st_centroid(sf::st_geometry(edges_sf)[path_labels$mid_edge])
)
edges_long <- path_df |>
  dplyr::mutate(path_id = dplyr::row_number()) |>
  tidyr::unnest(edge_paths) |>
  dplyr::mutate(
    from_name = locs$location_name[match(from, locs$network_nodes)],
    to_name = locs$location_name[match(to, locs$network_nodes)]
  )

edges_long_sf <- sf::st_sf(
  edges_long,
  geometry = sf::st_geometry(edges_sf)[edges_long$edge_paths]
)

sites <- sort(unique(edges_long_sf$from_name))

farms <- dplyr::filter(locs, location_type == "farm")
farm_box <- farms |>
  sf::st_geometry() |>
  sf::st_buffer(3000) |>
  sf::st_bbox()

box_poly <- sf::st_as_sfc(farm_box)

# plot all paths
all_paths_p <- ggplot() +
  # background
  geom_sf(data = bc_cropped) +
  #all paths
  geom_sf(
    data = network |>
      sfnetworks::activate("edges") |>
      dplyr::slice(unique(unlist(path_df$edge_paths))) |>
      sf::st_as_sf(),
    colour = "purple",
    alpha = 0.6,
    linewidth = 1
  ) +
  ggrepel::geom_label_repel(
    data = path_labels,
    aes(geometry = geometry, label = label),
    stat = "sf_coordinates",
    size = 3,
    box.padding = 0.6,
    min.segment.length = 0
  ) +
  geom_point(
    data = locs,
    aes(
      x = easting_x,
      y = northing_y,
      fill = location_type
    ),
    colour = "black",
    shape = 21,
    size = 4
  ) +
  theme_better()
ggsave(
  here::here("./paths-from-sites-to-farms.pdf"),
  all_paths_p
)

# 4.1 [PLOT PATHS DIFFERENT WAY] -----------------------------------------------

main_box <- sf::st_bbox(bc_cropped)
w <- as.numeric(main_box["xmax"] - main_box["xmin"])
h <- as.numeric(main_box["ymax"] - main_box["ymin"])

ins <- list(
  xmin = main_box["xmin"] + 0.63 * w,
  xmax = main_box["xmin"] + 0.98 * w,
  ymin = main_box["ymin"] + 0.51 * h,
  ymax = main_box["ymin"] + 0.98 * h
)

row_h <- 0.040 * h
tbl_x <- main_box["xmin"] + 0.045 * w
tbl_y0 <- main_box["ymin"] + 0.055 * h

corner_rows <- path_labels |>
  sf::st_drop_geometry() |>
  dplyr::group_by(from_name) |>
  dplyr::arrange(to_name, .by_group = TRUE) |>
  dplyr::mutate(
    row = dplyr::row_number(),
    txt = sprintf("%s: %.1f km", to_name, path_length / 1000)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(x = tbl_x, y = tbl_y0 + (max(row) - row) * row_h)

make_panel <- function(site) {
  e <- dplyr::filter(edges_long_sf, from_name == site)
  tb <- dplyr::filter(corner_rows, from_name == site)

  pts <- locs |>
    dplyr::mutate(
      plot_type = dplyr::if_else(
        location_type == "sampling" & location_name == site,
        "focal",
        location_type
      )
    )
  farms_ins <- farms |>
    dplyr::mutate(
      pt_fill = farm_cols[location_name] |>
        colorspace::lighten(0.4)
    )

  inset <- ggplot() +
    geom_sf(data = bc_cropped) +
    geom_sf(
      data = e,
      aes(colour = to_name),
      alpha = 0.9,
      linewidth = 2.4,
      lineend = "round",
      linejoin = "round"
    ) +
    geom_sf(
      data = farms_ins,
      aes(fill = I(pt_fill)),
      colour = "black",
      shape = 21,
      size = 6
    ) +
    ggrepel::geom_text_repel(
      data = farms,
      aes(geometry = geometry, label = location_name),
      stat = "sf_coordinates",
      size = 5,
      fontface = "bold",
      box.padding = 0.6,
      point.padding = 0.4,
      min.segment.length = 0,
      segment.size = 0.3,
      bg.colour = "white",
      bg.r = 0.15,
      seed = 1
    ) +
    scale_colour_manual(values = farm_cols, guide = "none") +
    coord_sf(
      xlim = c(farm_box["xmin"], farm_box["xmax"]),
      ylim = c(farm_box["ymin"], farm_box["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.8
      ),
      plot.margin = margin(2, 2, 2, 2)
    )

  ggplot() +
    geom_sf(data = bc_cropped) +
    geom_sf(
      data = e,
      aes(colour = to_name),
      alpha = 0.8,
      linewidth = 1,
      lineend = "round",
      linejoin = "round"
    ) +
    geom_sf(data = box_poly, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(
      data = pts,
      aes(fill = plot_type, size = plot_type),
      colour = "black",
      shape = 21
    ) +
    scale_size_manual(
      values = c("farm" = 4, "sampling" = 4, "focal" = 6),
      guide = "none"
    ) +
    annotation_custom(
      ggplotGrob(inset),
      xmin = ins$xmin,
      xmax = ins$xmax,
      ymin = ins$ymin,
      ymax = ins$ymax
    ) +
    annotate(
      "rect",
      xmin = tbl_x - 0.02 * w,
      xmax = tbl_x + 0.28 * w,
      ymin = tbl_y0 - 0.6 * row_h,
      ymax = tbl_y0 + (nrow(tb) - 1) * row_h + 0.6 * row_h,
      fill = "white",
      colour = "black",
      linewidth = 0.6
    ) +
    geom_text(
      data = tb,
      aes(x = x, y = y, label = txt, colour = to_name),
      hjust = 0,
      size = 4,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = farm_cols, name = "Farm") +
    scale_fill_manual(values = type_cols, name = "Location type") +
    coord_sf(expand = FALSE) +
    ggtitle(site) +
    theme_better()
}

faceted_paths <- wrap_plots(
  lapply(sites, make_panel),
  ncol = 3,
  guides = "collect"
)
faceted_paths_vert <- wrap_plots(
  lapply(sites, make_panel),
  ncol = 1,
  guides = "collect"
)

ggsave(
  here::here("paths-from-sites-to-farms-faceted.pdf"),
  faceted_paths,
  width = 26,
  height = 12
)
ggsave(
  here::here("paths-from-sites-to-farms-faceted.png"),
  faceted_paths_vert,
  width = 12,
  height = 30,
  dpi = 300
)
