# ============================================================
# FUNCTIONS FOR THE TIST LANDCOVER / ELIGIBILITY / FOREST-
# DISTANCE SECTION OF THE ERGM PIPELINE (_targets.R)
#
# Scope: everything the plan calls EXCEPT the network/ERGM
# functions (build_site_networks, fit_site_model,
# evaluate_overlap) -- those already exist elsewhere in your R/
# folder and are untouched, since that section is expensive and
# already has cached, verified results.
#
# This supersedes the earlier pipeline_functions.R for this
# plan. It corresponds to Sections A-E of the fuller v2 library
# built for the standalone analysis script, trimmed to only what
# this specific _targets.R plan actually calls, plus the small
# fidelity upgrades (colour palette, bracket-tag figure
# combination, more informative error message) pulled from the
# tested script.
# ============================================================

library(sf)
library(terra)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)


# ============================================================
# DATA LOADING
# ============================================================

load_TIST_data <- function(path) {
  
  if (!file.exists(path)) {
    stop("TIST data file not found: ", path)
  }
  
  tist_data <- utils::read.csv(path, stringsAsFactors = FALSE)
  
  if (nrow(tist_data) == 0) {
    stop("TIST data file contains no rows: ", path)
  }
  
  tist_data
}


load_uganda_boundary <- function(path) {
  
  sf::st_read(path, quiet = TRUE)
}


# ============================================================
# LAND-COVER MOSAIC / LEGEND / EXTRACTION / SUITABILITY
# ============================================================

create_landcover_legend <- function() {
  
  tibble::tibble(
    landcover_class = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
    landcover_name = c(
      "Tree cover", "Shrubland", "Grassland", "Cropland", "Built-up",
      "Bare / sparse vegetation", "Snow / ice", "Water bodies", "Wetlands"
    )
  )
}


build_landcover_mosaic_file <- function(files, boundary, output) {
  
  boundary <- terra::vect(boundary)
  boundary <- terra::project(boundary, "EPSG:4326")
  
  cropped <- lapply(files, function(f) {
    message("Processing: ", basename(f))
    r <- terra::rast(f)
    r <- terra::project(r, "EPSG:4326")
    terra::crop(r, boundary)
  })
  
  landcover <- do.call(terra::mosaic, c(cropped, fun = "first"))
  
  terra::writeRaster(landcover, output, overwrite = TRUE)
  
  # Returns a file path (not a SpatRaster) -- required for
  # tar_target(..., format = "file")
  output
}


extract_landcover <- function(tist_data, raster_file) {
  
  raster <- terra::rast(raster_file)
  
  pts <- terra::vect(
    tist_data,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  
  if (!terra::same.crs(raster, pts)) {
    pts <- terra::project(pts, raster)
  }
  
  extracted <- terra::extract(raster, pts)
  tist_data$landcover_class <- extracted[[2]]
  
  if (any(is.na(tist_data$landcover_class))) {
    warning(
      sum(is.na(tist_data$landcover_class)),
      " locations did not overlap landcover raster"
    )
  }
  
  tist_data
}


summarise_landcover <- function(tist_data, area, landcover_legend) {
  
  df <- tist_data |> dplyr::filter(Proj_Area == area)
  
  df |>
    dplyr::group_by(landcover_class) |>
    dplyr::summarise(
      n_groves = dplyr::n(),
      pct_groves = 100 * dplyr::n() / nrow(df),
      .groups = "drop"
    ) |>
    dplyr::left_join(landcover_legend, by = "landcover_class") |>
    dplyr::relocate(landcover_name, .after = landcover_class) |>
    dplyr::arrange(dplyr::desc(n_groves))
}


create_landcover_typology <- function() {
  
  tibble::tibble(
    landcover_name = c(
      "Tree cover", "Cropland", "Grassland", "Shrubland",
      "Built-up", "Bare / sparse vegetation"
    ),
    suitability_class = c(
      "Existing tree systems", "Open agricultural systems",
      "Open agricultural systems", "Semi-natural / extensifiable systems",
      "Excluded", "Excluded"
    ),
    eligible = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  )
}


# [Upgraded from earlier version] More informative unmapped-class
# error message, matching the tested standalone script.
apply_tree_suitability <- function(landcover_df, typology) {
  
  out <- landcover_df |> dplyr::left_join(typology, by = "landcover_name")
  
  if (any(is.na(out$suitability_class))) {
    stop(
      "Unmapped landcover classes detected: ",
      paste(unique(out$landcover_name[is.na(out$suitability_class)]), collapse = ", ")
    )
  }
  
  out
}


# ============================================================
# ELIGIBLE LAND / FOREST-RESERVE EXCLUSION
# (in-memory core logic + "_file" persistence wrappers)
# ============================================================

create_landcover_mask <- function(landcover_raster_file, eligible_classes = c(10, 20, 30, 40)) {
  
  r <- terra::rast(landcover_raster_file)
  
  if (terra::crs(r) == "") {
    stop("Land-cover raster has no CRS.")
  }
  
  eligible_mask <- r %in% eligible_classes
  eligible_mask[eligible_mask == 0] <- NA
  eligible_mask
}


create_landcover_mask_file <- function(
    landcover_raster_file,
    eligible_classes = c(10, 20, 30, 40),
    output = "Data/derived/eligible_landcover_mask.tif"
) {
  
  eligible_mask <- create_landcover_mask(landcover_raster_file, eligible_classes)
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(eligible_mask, output, overwrite = TRUE)
  output
}


prepare_forest_reserves_for_exclusion <- function(forest_reserves, landcover_raster_file) {
  
  r <- terra::rast(landcover_raster_file)
  raster_crs <- terra::crs(r, proj = TRUE)
  
  if (is.null(raster_crs) || is.na(raster_crs) || raster_crs == "") {
    stop("Land-cover raster has no valid CRS.")
  }
  
  forest_v <- terra::vect(forest_reserves)
  forest_crs <- terra::crs(forest_v, proj = TRUE)
  
  if (is.null(forest_crs) || is.na(forest_crs) || forest_crs == "") {
    stop("Forest-reserve layer has no valid CRS.")
  }
  
  forest_v <- terra::project(forest_v, raster_crs)
  forest_v <- terra::makeValid(forest_v)
  forest_v <- forest_v[!terra::is.empty(forest_v), ]
  
  if (nrow(forest_v) == 0) {
    stop("Forest-reserve layer contains no valid geometries.")
  }
  
  forest_v
}


prepare_forest_reserves_for_exclusion_file <- function(
    forest_reserves,
    landcover_raster_file,
    output = "Data/derived/forest_reserves_vector.gpkg"
) {
  
  forest_v <- prepare_forest_reserves_for_exclusion(forest_reserves, landcover_raster_file)
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  terra::writeVector(forest_v, output, overwrite = TRUE)
  output
}


# Forest reserves projected to a UTM CRS for distance work (used
# internally by prepare_forest_reserves_file()).
prepare_forest_distance_layer <- function(forest_reserves, utm_crs = 32636) {
  
  if (!inherits(forest_reserves, "sf")) {
    stop("forest_reserves must be an sf object")
  }
  
  forest_reserves |> sf::st_transform(utm_crs)
}


prepare_forest_reserves_file <- function(
    forest_reserves,
    target_crs = 32636,
    output = "Data/derived/forest_reserves_projected.gpkg"
) {
  
  forest_proj <- prepare_forest_distance_layer(forest_reserves, utm_crs = target_crs)
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(forest_proj, output, delete_dsn = TRUE, quiet = TRUE)
  output
}


create_forest_exclusion_mask <- function(forest_vector, landcover_raster_file) {
  
  r <- terra::rast(landcover_raster_file)
  
  if (terra::crs(r) == "") {
    stop("Land-cover raster has no CRS.")
  }
  
  if (!inherits(forest_vector, "SpatVector")) {
    forest_vector <- terra::vect(forest_vector)
  }
  
  if (terra::crs(forest_vector) == "") {
    stop("Forest-reserve vector has no CRS.")
  }
  
  if (terra::crs(forest_vector) != terra::crs(r)) {
    forest_vector <- terra::project(forest_vector, terra::crs(r))
  }
  
  forest_vector <- terra::makeValid(forest_vector)
  forest_vector <- forest_vector[!terra::is.empty(forest_vector), ]
  
  if (nrow(forest_vector) == 0) {
    stop("No valid forest-reserve geometries available for rasterisation.")
  }
  
  terra::rasterize(forest_vector, r, field = 1, background = NA)
}


create_forest_exclusion_mask_file <- function(
    forest_reserves_vector_file,
    landcover_raster_file,
    output = "Data/derived/forest_mask.tif"
) {
  
  forest_v <- terra::vect(forest_reserves_vector_file)
  reserve_mask <- create_forest_exclusion_mask(forest_v, landcover_raster_file)
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(reserve_mask, output, overwrite = TRUE)
  output
}


exclude_forest_reserves <- function(landcover_mask, forest_mask) {
  
  landcover_mask <- terra::rast(landcover_mask)
  forest_mask <- terra::rast(forest_mask)
  
  if (!terra::compareGeom(landcover_mask, forest_mask, stopOnError = FALSE)) {
    stop("Landcover mask and forest mask do not have matching geometry.")
  }
  
  terra::mask(landcover_mask, forest_mask, inverse = TRUE)
}


exclude_forest_reserves_file <- function(
    eligible_mask_file,
    forest_mask_file,
    output = "Data/derived/eligible_area.tif"
) {
  
  eligible <- exclude_forest_reserves(eligible_mask_file, forest_mask_file)
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(eligible, output, overwrite = TRUE)
  output
}


# ============================================================
# PROJECT AREAS
# ============================================================

create_project_area <- function(districts, district_names) {
  
  area <- districts |> dplyr::filter(ADM2_EN %in% district_names)
  
  if (nrow(area) == 0) {
    stop("No matching districts found.")
  }
  
  sf::st_as_sf(sf::st_union(area))
}


# Tolerant of eligible_raster being either a file path (string)
# or an already-loaded SpatRaster -- lets this be called with
# eligible_area_file directly.
mask_project_area <- function(eligible_raster, project_area) {
  
  if (is.character(eligible_raster)) {
    eligible_raster <- terra::rast(eligible_raster)
  }
  
  area_v <- terra::vect(project_area)
  area_v <- terra::project(area_v, terra::crs(eligible_raster))
  area_v <- terra::makeValid(area_v)
  
  eligible_raster |>
    terra::crop(area_v) |>
    terra::mask(area_v)
}


# [NEW] File-persisting wrapper. SpatRaster objects hold an
# external C++ pointer -- storing one directly as a target value
# (no format = "file") works only within a single live R session;
# as soon as targets serialises it to disk and a later target
# reloads it (a new process, or even the same session after a
# round trip through storage), that pointer is dead and any
# terra:: call on it fails with "external pointer is not valid".
# This wrapper writes to disk immediately and returns the path,
# matching every other raster target in this pipeline
# (landcover_raster_file, eligible_mask_file, forest_mask_file,
# eligible_area_file, ...). mask_project_area() itself is
# unchanged and still useful for interactive/in-session use.
mask_project_area_file <- function(eligible_raster, project_area, output) {
  
  masked <- mask_project_area(eligible_raster, project_area)
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(masked, output, overwrite = TRUE)
  output
}


# ============================================================
# TIST / RANDOM POINTS AND FOREST-DISTANCE ANALYSIS
# ============================================================

prepare_TIST_distance_points <- function(tist_data, area, utm_crs = 32636) {
  
  points <- tist_data |> dplyr::filter(Proj_Area == area)
  
  if (nrow(points) == 0) {
    stop("No TIST locations found for: ", area)
  }
  
  points |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(utm_crs)
}


# Tolerant of eligible_raster being a file path or SpatRaster.
# [FIXED] Now reprojects to utm_crs before returning, matching
# prepare_TIST_distance_points(). Without this, points came back
# in the landcover raster's CRS (EPSG:4326) while
# forest_reserves_projected_file is in UTM 36N (EPSG:32636) --
# nearest_forest_distance()'s st_distance() call then fails with
# "st_crs(x) == st_crs(y) is not TRUE".
prepare_random_ag_points <- function(eligible_raster, n_points = 1000, seed = 123, utm_crs = 32636) {
  
  if (is.character(eligible_raster)) {
    eligible_raster <- terra::rast(eligible_raster)
  }
  
  set.seed(seed)
  
  points <- terra::spatSample(
    eligible_raster, size = n_points,
    method = "random", as.points = TRUE, na.rm = TRUE
  )
  
  sf::st_as_sf(points) |> sf::st_transform(utm_crs)
}


# Tolerant of forest_layer being a file path (string) or an
# already-loaded sf object -- lets this be called with
# forest_reserves_projected_file directly.
# [DEFENSIVE FIX] Reprojects points to forest_layer's CRS if they
# don't already match, rather than letting st_distance() error --
# guards against any future caller of this function forgetting to
# pre-transform its points (as prepare_random_ag_points() did).
nearest_forest_distance <- function(points, forest_layer) {
  
  if (is.character(forest_layer)) {
    forest_layer <- sf::st_read(forest_layer, quiet = TRUE)
  }
  
  if (nrow(forest_layer) == 0) {
    stop("Forest reserve layer contains no features")
  }
  
  if (!identical(sf::st_crs(points), sf::st_crs(forest_layer))) {
    points <- sf::st_transform(points, sf::st_crs(forest_layer))
  }
  
  dist_matrix <- sf::st_distance(points, forest_layer)
  points$Dist_To_Forest_m <- apply(dist_matrix, 1, min) |> as.numeric()
  points
}


create_distance_comparison <- function(observed, random) {
  
  observed_df <- observed |>
    sf::st_drop_geometry() |>
    dplyr::select(Dist_To_Forest_m) |>
    dplyr::mutate(type = "TIST")
  
  random_df <- random |>
    sf::st_drop_geometry() |>
    dplyr::select(Dist_To_Forest_m) |>
    dplyr::mutate(type = "Random_ag")
  
  dplyr::bind_rows(observed_df, random_df) |>
    dplyr::filter(!is.na(Dist_To_Forest_m), Dist_To_Forest_m >= 0)
}


summarise_forest_distance <- function(comparison) {
  
  comparison |>
    dplyr::group_by(type) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(Dist_To_Forest_m, na.rm = TRUE),
      median_distance = median(Dist_To_Forest_m, na.rm = TRUE),
      q25 = quantile(Dist_To_Forest_m, 0.25, na.rm = TRUE),
      q75 = quantile(Dist_To_Forest_m, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}


test_forest_distance <- function(comparison) {
  
  wilcox.test(Dist_To_Forest_m ~ type, data = comparison, alternative = "less")
}


# ============================================================
# FIGURES
# ============================================================

# The colour palette verified in the tested standalone script.
get_landcover_colors <- function() {
  
  c(
    "Tree cover" = "#1b7837",
    "Shrubland" = "#7fbf7b",
    "Grassland" = "#d9f0d3",
    "Cropland" = "#dfc27d",
    "Built-up" = "#b2182b",
    "Bare / sparse vegetation" = "#bababa"
  )
}


# [Upgraded from earlier version] Accepts an optional fill_cols
# palette (e.g. get_landcover_colors()); falls back to ggplot's
# default palette if omitted, so existing 2-argument calls still
# work unchanged.
plot_landcover_pie <- function(df, area_name, fill_cols = NULL) {
  
  df <- df |> dplyr::mutate(pct_label = round(pct_groves, 1))
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = "", y = pct_groves, fill = landcover_name)
  ) +
    ggplot2::geom_col(width = 1, colour = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(pct_label, "%")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.5
    ) +
    ggplot2::labs(
      title = paste0(area_name, ": TIST grove land cover"),
      fill = "Land cover"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(fill_cols)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_cols)
  }
  
  p
}


plot_suitability_bar <- function(df, area_name) {
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = reorder(suitability_class, pct_groves), y = pct_groves)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = "Percentage of groves",
      title = paste0(area_name, ": planting suitability classes")
    ) +
    ggplot2::theme_minimal()
}


# [Upgraded from earlier version] Bracketed sequential tags --
# "(a)", "(b)" -- matching the tested standalone script's
# combine_landcover_piecharts(), kept under the name this plan
# actually calls (combine_two_area_plots) so no target needs
# renaming.
combine_two_area_plots <- function(plot1, plot2) {
  
  (plot1 + plot2 +
     patchwork::plot_layout(ncol = 2) +
     patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 16, face = "bold"),
      plot.tag.position = "top"
    )
}
