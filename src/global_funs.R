##' File Description
##' AUTHOR: Cole B. Brookson
##' DATE OF CREATION: 2022-10-14
#'
#' This targets file contains all functions that don't relate to a specific
#' part of the analysis, but are required to perform general tasks
#'
#' All functions are documented using the roxygen2 framework and the docstring
#' library
#'

`%notin%` <- Negate(`%in%`)

#' Foundation Theme
#'
#' This theme is designed to be a foundation from which to build new
#' themes, and not meant to be used directly. \code{theme_foundation()}
#' is a complete theme with only minimal number of elements defined.
#' It is easier to create new themes by extending this one rather
#' than \code{\link[ggplot2]{theme_gray}()} or \code{\link[ggplot2]{theme_bw}()},
#' because those themes define elements deep in the hierarchy.
#'
#' This theme takes \code{\link[ggplot2]{theme_gray}()} and sets all
#' \code{colour} and \code{fill} values to \code{NULL}, except for the top-level
#' elements (\code{line}, \code{rect}, and \code{title}), which have
#' \code{colour = "black"}, and \code{fill = "white"}. This leaves the spacing
#' and-non colour defaults of the default \pkg{ggplot2} themes in place.
#'
#' @export
#' @inheritParams ggplot2::theme_grey
#'
#' @family themes
#' @export
#' @importFrom ggplot2 theme_grey
theme_foundation <- function(base_size = 16, base_family = "") {
    thm <- theme_grey(base_size = base_size, base_family = base_family)
    for (i in names(thm)) {
        if ("colour" %in% names(thm[[i]])) {
            thm[[i]]["colour"] <- list(NULL)
        }
        if ("fill" %in% names(thm[[i]])) {
            thm[[i]]["fill"] <- list(NULL)
        }
    }
    thm + theme(
        panel.border = element_rect(fill = NA),
        legend.background = element_rect(colour = NA),
        line = element_line(colour = "black"),
        rect = element_rect(fill = "white", colour = "black"),
        text = element_text(colour = "black")
    )
}

#' Theme Base
#'
#' Theme similar to the default settings of the \sQuote{base} R graphics.
#'
#' @inheritParams ggplot2::theme_bw
#' @export
#' @family themes
#' @example inst/examples/ex-theme_base.R
theme_base <- function(base_size = 16, base_family = "") {
    theme_foundation() +
        theme(
            line = element_line(
                colour = "black",
                lineend = "round",
                linetype = "solid"
            ),
            rect = element_rect(
                fill = "white",
                colour = "black",
                linetype = "solid"
            ),
            text = element_text(
                colour = "black",
                face = "plain",
                family = base_family,
                size = base_size,
                vjust = 0.5,
                hjust = 0.5,
                lineheight = 1
            ),
            panel.grid = element_blank(),
            strip.background = element_rect(colour = NA),
            legend.key = element_rect(colour = NA),
            title = element_text(size = rel(1)),
            plot.title = element_text(size = rel(1.2), face = "bold"),
            strip.text = element_text(),
            axis.ticks.length = unit(0.5, "lines"),
            # add my addition here
            plot.background = element_rect(colour = NA)
        )
    # TODO: get margins right
}

#' Compute total length of a path from edge indices
#'
#' @description
#' Given a set of edge indices from a network path, this function computes the
#' total length of the corresponding path in the spatial network. The result is
#' returned as a numeric value with units dropped.
#'
#' @param net An `sfnetwork` object representing the spatial network.
#' @param temp_edges A vector of edge indices representing a path through the network.
#'
#' @return A numeric scalar giving the total length of the path (in meters).
#' @export
slice_fun <- function(net, temp_edges) {
    return(net |>
        tidygraph::activate("edges") |>
        dplyr::slice(temp_edges) |>
        sf::st_as_sf() |>
        sf::st_combine() |>
        sf::st_length() |>
        units::drop_units())
}

#' Compute all shortest paths from one node to a set of target nodes
#'
#' @description
#' Uses `sfnetworks::st_network_paths()` to compute shortest paths from a single
#' node in a spatial network to a set of target nodes, with optional weighting.
#'
#' @param from_node An integer index of the source node in the network.
#' @param to_nodes A vector of integer indices representing target nodes.
#' @param net An `sfnetwork` object representing the spatial network.
#'
#' @return A list object returned by `sfnetworks::st_network_paths()` that includes
#'         node and edge paths.
#' @export
get_paths_from_one <- function(from_node, to_nodes, net) {
    sfnetworks::st_network_paths(
        x = net,
        from = from_node,
        to = to_nodes,
        weights = "weight"
    )
}

#' Compute pairwise shortest path distances between nodes in an sfnetwork
#'
#' @description
#' Given a set of node indices from an `sfnetwork`, this function computes the
#' shortest path distances between all unique pairs using `sfnetworks::st_network_paths()`
#' and a helper to compute total path length. Supports optional parallel execution with `furrr`.
#'
#' @param node_ids An integer vector of node indices corresponding to points of interest.
#' @param net An `sfnetwork` object representing the spatial network with weighted edges.
#' @param parallel Logical. If `TRUE`, uses `furrr::future_map()` for parallel computation.
#'
#' @return A `tibble` with columns:
#'   - `from`: source node ID
#'   - `to`: target node ID
#'   - `edge_paths`: list of edge indices
#'   - `path_length`: total length of path
#' @export
get_pairwise_network_distances <- function(
    node_ids,
    net,
    parallel = TRUE) {
    # choose map function based on parallel flag
    compute_paths <- if (parallel) furrr::future_map else purrr::map

    # compute all shortest paths
    all_paths <- compute_paths(node_ids, function(from_node) {
        get_paths_from_one(from_node, node_ids, net)
    },
    .progress = TRUE
    )

    # flatten into dataframe
    path_df <- purrr::imap_dfr(
        all_paths,
        ~ tibble::tibble(
            from = node_ids[.y],
            to = node_ids,
            edge_paths = .x$edge_paths
        ),
        .progress = TRUE
    )

    # compute path lengths
    path_df <- path_df |>
        dplyr::mutate(
            path_length = furrr::future_map_dbl(
                edge_paths, ~ slice_fun(net, .x),
                seed = TRUE,
                .progress = TRUE
            )
        ) |>
        dplyr::filter(from != to)

    return(path_df)
}


# #' Compute pairwise shortest path distances between nodes in an sfnetwork
# #'
# #' @description
# #' Given a set of node indices from an `sfnetwork`, this function computes the
# #' shortest path distances between all unique pairs of nodes using
# #' `sfnetworks::st_network_paths()` and a helper function to sum the path lengths.
# #' A progress bar is included for long-running operations.
# #'
# #' @param node_ids An integer vector of node indices corresponding to points of interest
# #'        (e.g., nearest nodes to sample locations).
# #' @param net An `sfnetwork` object representing the spatial network, with a
# #'        `weight` column on the edges (e.g., length in meters).
# #'
# #' @return A `tibble` with columns:
# #'   - `from`: source node ID
# #'   - `to`: target node ID
# #'   - `path_length`: numeric shortest path length between `from` and `to`
# #' @export
# get_pairwise_network_distances <- function(node_ids, net) {
#     # # setup progress handler
#     # progressr::handlers(global = TRUE)
#     # p <- progressr::progressor(along = node_ids)

#     all_paths <- purrr::map(
#         node_ids,
#         function(from_node) {
#             print(sprintf("Calculating paths from node %s", from_node))
#             get_paths_from_one(from_node, node_ids, net)
#         },
#         .progress = list(
#             clear = FALSE,
#             format_done = "1/3 tasks completed",
#             name = "Getting All Paths (1/3 tasks)"
#         )
#     )

#     # Flatten into a dataframe
#     path_df <- purrr::imap_dfr(
#         all_paths,
#         ~ tibble::tibble(
#             from = node_ids[.y],
#             to = node_ids,
#             edge_paths = .x$edge_paths
#         )
#     )

#     # Compute lengths using slice_fun
#     path_df <- path_df |>
#         dplyr::mutate(
#             path_length = purrr::map_dbl(edge_paths, ~ slice_fun(net, .x))
#         ) |>
#         dplyr::filter(from != to)

#     return(path_df)
# }

#' Plot paths from a given node in a spatial network with labeling and optional title
#'
#' @description
#' Plots spatial paths from a single source node to:
#' - all other nodes (`to_nodes = "all"`),
#' - the two nearest nodes (`to_nodes = "nearest"`),
#' - or a specific target node or vector of targets (`to_nodes = <node ID(s)>`).
#' Adds node ID labels and an optional title for distance.
#'
#' @param from_node Integer. The node ID to plot paths **from**.
#' @param to_nodes Either:
#'   - a single integer or vector of node IDs,
#'   - the string `"all"` to plot all paths from `from_node`,
#'   - or the string `"nearest"` to plot paths to the 2 closest nodes.
#' @param network An `sfnetwork` object representing the spatial network.
#' @param distance_table A data.frame or tibble with columns `from`, `to`, `path_length`, and `edge_paths`.
#' @param background_sf An `sf` object representing the map background (e.g., land or water mask).
#' @param samples_sf An `sf` object of sampling locations, must include a `nearest_node` column.
#' @param farms_sf An `sf` object of farm locations.
#' @param zoom_to_extent Logical. If `TRUE`, zooms the plot extent to just the region around the plotted paths (+20% padding).
#'
#' @return A `ggplot` object showing the selected paths and labeled points.
#' @export
plot_paths_from_node <- function(
    from_node, to_nodes, network, distance_table, background_sf,
    samples_sf, farms_sf, zoom_to_extent = FALSE) {
    # validate `to_nodes`
    if (!(is.numeric(to_nodes) || to_nodes %in% c("all", "nearest"))) {
        stop('`to_nodes` must be numeric, "all", or "nearest"')
    }

    # filter paths from the specified source node
    paths_from <- dplyr::filter(distance_table, from == from_node)

    # handle selection
    if (is.numeric(to_nodes)) {
        paths_plot <- dplyr::filter(paths_from, to %in% to_nodes)
    } else if (to_nodes == "all") {
        paths_plot <- paths_from
    } else if (to_nodes == "nearest") {
        paths_plot <- paths_from |>
            dplyr::mutate(length = sf::st_length(geometry)) |>
            dplyr::arrange(length) |>
            dplyr::slice(1:2)
    }

    # label points
    sample_labels <- samples_sf |>
        dplyr::mutate(label = as.character(nearest_node)) |>
        dplyr::filter(label %in% c(paths_plot$from, paths_plot$to))

    # extract edge geometries
    edges_sf <- network |>
        sfnetworks::activate("edges") |>
        dplyr::slice(unlist(paths_plot$edge_paths)) |>
        sf::st_as_sf()

    # base plot
    p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = background_sf) +
        ggplot2::geom_sf(data = edges_sf, color = "orange", size = 0.6) +
        ggplot2::geom_sf(
            data = samples_sf |>
                dplyr::filter(nearest_node %in% c(paths_plot$from, paths_plot$to)),
            color = "blue", size = 5
        ) +
        ggplot2::geom_sf(
            data = farms_sf, shape = 21, fill = "red",
            color = "black", size = 1.5
        ) +
        ggplot2::geom_sf_text(
            data = sample_labels,
            ggplot2::aes(label = label), size = 4, nudge_y = 200
        ) +
        theme_base()

    # optional: add title if only one distance
    if (is.numeric(to_nodes) && length(to_nodes) == 1) {
        path_dist <- min(paths_plot$path_length)
        p <- p + ggplot2::ggtitle(sprintf("Distance: %.1f meters", path_dist))
    }

    # optional: zoom to region of interest
    if (zoom_to_extent) {
        # get extent of edge geometries
        bbox <- sf::st_bbox(edges_sf)
        x_pad <- 0.2 * (bbox["xmax"] - bbox["xmin"])
        y_pad <- 0.2 * (bbox["ymax"] - bbox["ymin"])

        p <- p +
            ggplot2::coord_sf(
                xlim = c(bbox["xmin"] - x_pad, bbox["xmax"] + x_pad),
                ylim = c(bbox["ymin"] - y_pad, bbox["ymax"] + y_pad),
                expand = FALSE
            )
    }

    return(p)
}
