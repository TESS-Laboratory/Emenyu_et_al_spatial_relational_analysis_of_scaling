################################################################################
# TIST Uganda — Land Cover Overlay, Planting-Suitability Sampling Frame,
# and Forest-Proximity Comparison (Observed Groves vs Random Agricultural Land)
################################################################################

## ---- Packages --------------------------------------------------------------
## Trimmed to the packages whose functions are actually called below:
##   terra      - raster/vector handling, mosaicking, extraction, masking, sampling
##   sf         - vector geometry, projections, distance calculations
##   tidyverse  - dplyr/ggplot2/tibble/purrr workflow
##   patchwork  - plot_layout() / plot_annotation() / "+" combining of ggplots
##   magick     - final image stacking/annotation of composite figures
library(terra)
library(sf)
library(tidyverse)
library(patchwork)
library(magick)


################################################################################
# 1. LOAD DATA & IDENTIFY PROJECT-AREA DISTRICTS
################################################################################

TistDist_Dat <- read.csv("Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity.csv") %>%
  drop_na() %>%
  distinct()

Bushenyi_districts <- TistDist_Dat %>%
  filter(Proj_Area == "Bushenyi") %>%
  distinct(Admin_Districts) %>%
  arrange(Admin_Districts)
Bushenyi_districts

Soroti_districts <- TistDist_Dat %>%
  filter(Proj_Area == "Soroti") %>%
  distinct(Admin_Districts) %>%
  arrange(Admin_Districts)
Soroti_districts


################################################################################
# 2. DISTRICT & FOREST RESERVE BOUNDARIES
################################################################################

districts_sf <- st_read("C:/Users/ae474/OneDrive - University of Exeter/PhD work in progress/Data/R Large Data/UGshapefiles/uga_admbnda_adm2_ubos_20200824.shp")

## Bushenyi project area districts
Bushenyi_names <- c(
  "Buhweju", "Bushenyi", "Ibanda", "Kamwenge", "Kasese",
  "Kitagwenda", "Kyenjojo", "Mitooma", "Ntungamo",
  "Rubirizi", "Rwampara", "Sheema"
)

## Soroti project area districts
Soroti_names <- c(
  "Alebtong", "Amuria", "Kalaki", "Kapelebyong", "Serere", "Soroti"
)

## Filter from data
Bushenyi_sf <- districts_sf |> filter(ADM2_EN %in% Bushenyi_names)
Soroti_sf   <- districts_sf |> filter(ADM2_EN %in% Soroti_names)

## Dissolve to project areas (still sf)
Bushenyi_area_sf <- st_as_sf(st_union(Bushenyi_sf))
Soroti_area_sf   <- st_as_sf(st_union(Soroti_sf))

forest_reserves <- st_read("Data/Uganda_Forest_Reserves/Uganda_Forest_Reserves.shp")


################################################################################
# 3. BUILD THE NATIONAL LAND COVER MOSAIC (ESA WorldCover 2020)
################################################################################

Uganda_sf <- st_read("C:/Users/ae474/OneDrive - University of Exeter/PhD work in progress/Data/R Large Data/UGshapefiles/uga_admbnda_adm0_ubos_20200824.shp")
Uganda_sf

## Convert to vector in terra
Uganda <- vect(Uganda_sf)
crs(Uganda)  # has to be EPSG:4326 for consistent distance comparisons later

Uganda <- project(Uganda, "EPSG:4326")
crs(Uganda)

files <- list.files(
  "Data/landcover/ESA_WorldCover_2020/",
  pattern = "_Map.tif$",
  full.names = TRUE
)

## Only keep files covering Uganda
uganda_tiles <- files[grepl("N0[03]|S03", files) & grepl("E03[03-6]", files)]
uganda_tiles

## Start with tiles that most definitely cover Uganda
uganda_tiles <- c(
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_N00E027_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_N00E030_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_N00E033_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_N03E027_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_N03E030_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_N03E033_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_S03E027_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_S03E030_Map.tif",
  "Data/landcover/ESA_WorldCover_2020/ESA_WorldCover_10m_2020_v100_S03E033_Map.tif"
)

cropped <- list()

for (f in uganda_tiles) {
  r <- rast(f)
  r <- project(r, "EPSG:4326")       # project raster to EPSG:4326 (safe)
  cropped[[length(cropped) + 1]] <- crop(r, Uganda)  # crop to Uganda
}

## Merge cropped tiles
landcover_Uganda <- do.call(mosaic, c(cropped, fun = "first"))

windows(width = 9, height = 9)  ## figure margins too large hence the use of windows
terra::plot(landcover_Uganda)
plot(Uganda, add = TRUE, border = "red", lwd = 2)

## Save to file
png("Output/Manuscript 3 graphs/uganda_landcover.png", width = 9, height = 9, units = "in", res = 300)
terra::plot(landcover_Uganda)
plot(Uganda, add = TRUE, border = "red", lwd = 2)
dev.off()

## Save mosaic for later use
writeRaster(
  landcover_Uganda,
  "Data/landcover/ESA_WorldCover_2020/Uganda_mosaic.tif",
  overwrite = TRUE
)


################################################################################
# 4. OVERLAY TIST GROVES ON LAND COVER & SUMMARISE LAND-USE CLASSES
################################################################################

landcover   <- landcover_Uganda  # SpatRaster
tist_points <- TistDist_Dat      # farmer locations

## Land cover legend table
landcover_legend <- tibble::tibble(
  landcover_class = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
  landcover_name = c(
    "Tree cover",
    "Shrubland",
    "Grassland",
    "Cropland",
    "Built-up",
    "Bare / sparse vegetation",
    "Snow / ice",
    "Water bodies",
    "Wetlands"
  )
)

## Convert points into a vector
tist_vect <- vect(tist_points, geom = c("longitude", "latitude"), crs = crs(landcover))

## Extract landcover value for each point
tist_points$landcover_class <- terra::extract(landcover, tist_vect)[, 2]

## Summarise number of groves per landcover class (national)
landcover_summary <- tist_points %>%
  group_by(landcover_class) %>%
  summarise(
    n_groves   = n(),
    pct_groves = 100 * n() / nrow(tist_points),
    .groups = "drop"
  ) %>%
  left_join(landcover_legend, by = "landcover_class") %>%
  relocate(landcover_name, .after = landcover_class)

landcover_summary

## ---- Helper: landcover summary for a given project area --------------------

landcover_summary_by_area <- function(data, raster, area_filter) {

  df <- data %>% filter({{ area_filter }})

  pts <- vect(df, geom = c("longitude", "latitude"), crs = crs(raster))

  df$landcover_class <- terra::extract(raster, pts)[, 2]

  df %>%
    group_by(landcover_class) %>%
    summarise(
      n_groves   = n(),
      pct_groves = 100 * n() / nrow(df),
      .groups = "drop"
    ) %>%
    left_join(landcover_legend, by = "landcover_class") %>%
    relocate(landcover_name, .after = landcover_class) %>%
    arrange(desc(n_groves))
}

bushenyi_landcover <- landcover_summary_by_area(
  data = TistDist_Dat, raster = landcover_Uganda, area_filter = Proj_Area == "Bushenyi"
)
bushenyi_landcover

soroti_landcover <- landcover_summary_by_area(
  data = TistDist_Dat, raster = landcover_Uganda, area_filter = Proj_Area == "Soroti"
)
soroti_landcover


################################################################################
# 5. LAND COVER PIE CHARTS
################################################################################

plot_landcover_pie <- function(df, district_name, fill_cols = NULL) {

  df <- df %>% mutate(pct_lab = round(pct_groves, 1))

  p <- ggplot(df, aes(x = "", y = pct_groves, fill = landcover_name)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = paste0(pct_lab, "%")),
      position = position_stack(vjust = 0.5),
      size = 3.5
    ) +
    labs(
      title = paste0(district_name, ": TIST grove Land Cover Distribution"),
      fill  = "Land cover"
    ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  if (!is.null(fill_cols)) {
    p <- p + scale_fill_manual(values = fill_cols)
  }

  return(p)
}

lc_cols <- c(
  "Tree cover"               = "#1b7837",
  "Shrubland"                = "#7fbf7b",
  "Grassland"                = "#d9f0d3",
  "Cropland"                 = "#dfc27d",
  "Built-up"                 = "#b2182b",
  "Bare / sparse vegetation" = "#bababa"
)

bushenyi_plot <- plot_landcover_pie(bushenyi_landcover, district_name = "Bushenyi", fill_cols = lc_cols)
bushenyi_plot

soroti_plot <- plot_landcover_pie(soroti_landcover, district_name = "Soroti", fill_cols = lc_cols)
soroti_plot

BushSrtgrove_Londcoverdtn_piechart <- bushenyi_plot + soroti_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "a",   # sequential letters
    tag_prefix = "(",   # add opening bracket
    tag_suffix = ")"    # add closing bracket
  ) &
  theme(
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = "top"  # put labels above the plots
  )

BushSrtgrove_Londcoverdtn_piechart

ggsave("Output/Manuscript 3 graphs/BushSrtgrove_Londcoverdtn_piechart.png",
       BushSrtgrove_Londcoverdtn_piechart, dpi = 600, width = 12, height = 6)


################################################################################
# 6. TREE-PLANTING SUITABILITY TYPOLOGY & SAMPLING-SPACE DEFINITION
################################################################################
#
# NOTE: In this work, land-use classes are treated as indicators of
# biophysical surface characteristics rather than direct land use.
#
# Based on the land-use classes to which observed grove locations were
# matched, we define a priority land-cover typology for tree-planting
# suitability.
################################################################################

landcover_typology <- tibble(
  landcover_name = c(
    "Tree cover",
    "Cropland",
    "Grassland",
    "Shrubland",
    "Built-up",
    "Bare / sparse vegetation"
  ),

  # Conceptual tree-planting suitability class
  suitability_class = c(
    "Existing tree systems",
    "Open agricultural systems",
    "Open agricultural systems",
    "Semi-natural / extensifiable systems",
    "Excluded",
    "Excluded"
  ),

  # Explicit eligibility flag for analytic sample space
  eligible = c(
    TRUE,   # Tree cover
    TRUE,   # Cropland
    TRUE,   # Grassland
    TRUE,   # Shrubland
    FALSE,  # Built-up
    FALSE   # Bare / sparse vegetation
  )
)

## ---- Apply tree-planting suitability typology to landcover summary data ----

apply_tree_suitability <- function(landcover_df) {

  out <- landcover_df %>%
    left_join(landcover_typology, by = "landcover_name")

  # Safety check: ensure all classes were mapped
  if (any(is.na(out$suitability_class))) {
    stop(
      "Unmapped landcover classes detected: ",
      paste(unique(out$landcover_name[is.na(out$suitability_class)]), collapse = ", ")
    )
  }

  return(out)
}

bushenyi_suitability <- apply_tree_suitability(bushenyi_landcover)
bushenyi_suitability

soroti_suitability <- apply_tree_suitability(soroti_landcover)
soroti_suitability

## ---- Share of groves excluded vs eligible -----------------------------------

bushenyi_suitability %>%
  group_by(eligible) %>%
  summarise(n_groves = sum(n_groves), pct_groves = sum(pct_groves))
## eligible groves: 2814 (99.1%); excluded: 25

soroti_suitability %>%
  group_by(eligible) %>%
  summarise(n_groves = sum(n_groves), pct_groves = sum(pct_groves))
## eligible groves: 2315 (98.3%); excluded: 41

## ---- Aggregate suitability classes -------------------------------------------

bushenyi_by_suitability <- bushenyi_suitability %>%
  filter(eligible) %>%    # explicit exclusion
  group_by(suitability_class) %>%
  summarise(n_groves = sum(n_groves), pct_groves = sum(pct_groves), .groups = "drop")
bushenyi_by_suitability

soroti_by_suitability <- soroti_suitability %>%
  filter(eligible) %>%
  group_by(suitability_class) %>%
  summarise(n_groves = sum(n_groves), pct_groves = sum(pct_groves), .groups = "drop")
soroti_by_suitability


################################################################################
# 7. BUILD ELIGIBILITY RASTER (LAND COVER + FOREST-RESERVE EXCLUSION)
################################################################################

eligible_lc_classes <- c(
  10,  # Tree cover
  20,  # Shrubland
  30,  # Grassland
  40   # Cropland
)

## Eligibility mask (TRUE for eligible classes, NA otherwise)
eligible_mask1 <- landcover_Uganda %in% eligible_lc_classes
eligible_mask1[eligible_mask1 == 0] <- NA  ## prevents non-eligible area from being a valid sampling area

## Forest reserves in raster CRS
forest_reserves_v <- vect(forest_reserves) |>
  project(crs(landcover_Uganda)) |>
  makeValid()

## Rasterize forest reserves as an exclusion mask
forest_reserve_mask <- rasterize(
  forest_reserves_v,
  landcover_Uganda,
  field      = 1,
  background = NA
)

## Apply exclusion to eligibility mask
eligible_mask <- mask(
  eligible_mask1,
  forest_reserve_mask,
  inverse = TRUE
)

## Mask eligible areas to Bushenyi and Soroti project areas
Bushenyi_area <- vect(Bushenyi_area_sf) |>
  project(crs(landcover_Uganda)) |>
  makeValid()

Soroti_area <- vect(Soroti_area_sf) |>
  project(crs(landcover_Uganda)) |>
  makeValid()

eligible_bushenyi <- mask(crop(eligible_mask, Bushenyi_area), Bushenyi_area)
eligible_soroti   <- mask(crop(eligible_mask, Soroti_area),   Soroti_area)

writeRaster(eligible_bushenyi, "Output/Manuscript 3 graphs/eligible_bushenyi.tif", overwrite = TRUE)
writeRaster(eligible_soroti,   "Output/Manuscript 3 graphs/eligible_soroti.tif",   overwrite = TRUE)

## ---- Sanity check plots -------------------------------------------------------

plot_eligibility <- function(r, area, border_col) {
  terra::plot(r)
  plot(area, add = TRUE, border = border_col, lwd = 2)
}

windows(width = 9, height = 9)
plot_eligibility(eligible_bushenyi, Bushenyi_area, "red")

windows(width = 9, height = 9)
plot_eligibility(eligible_soroti, Soroti_area, "blue")

png("Output/Manuscript 3 graphs/eligible_bushenyi.png", width = 9, height = 9, units = "in", res = 300)
plot_eligibility(eligible_bushenyi, Bushenyi_area, "red")
dev.off()

png("Output/Manuscript 3 graphs/eligible_soroti.png", width = 9, height = 9, units = "in", res = 300)
plot_eligibility(eligible_soroti, Soroti_area, "blue")
dev.off()


################################################################################
# 8. RANDOM SAMPLE OF ELIGIBLE AGRICULTURAL LAND
################################################################################
#
# Are groves preferentially located near forest reserves, or just positioned
# where agricultural land happens to be available?
#
# H0: There is no difference between the distance of groves to the nearest
#     forest reserve and the distance of random agricultural land to the
#     nearest forest reserve.
################################################################################

set.seed(123)   # reproducible
n_points <- 1000  # adjust as needed

points_bushenyi <- spatSample(
  eligible_bushenyi,
  size    = n_points,
  method  = "random",
  as.points = TRUE,
  na.rm   = TRUE
)

points_soroti <- spatSample(
  eligible_soroti,
  size    = n_points,
  method  = "random",
  as.points = TRUE,
  na.rm   = TRUE
)

## ---- Sanity check sampling plots ----------------------------------------------

plot_sampling <- function(points, area, outfile, border_col) {
  png(outfile, width = 9, height = 9, units = "in", res = 300)
  terra::plot(points)
  plot(area, add = TRUE, border = border_col, lwd = 2)
  dev.off()
}

plot_sampling(points_bushenyi, Bushenyi_area,
              "Output/Manuscript 3 graphs/samplingpoints_bushenyi.png", border_col = "red")

plot_sampling(points_soroti, Soroti_area,
              "Output/Manuscript 3 graphs/samplingpoints_soroti.png", border_col = "blue")

## Final validation
stopifnot(all(relate(points_bushenyi, Bushenyi_area, "within")))
stopifnot(all(relate(points_soroti,   Soroti_area,   "within")))

dev.off()  # turn off external windows plotting area


################################################################################
# 9. OBSERVED (TIST) vs RANDOM AGRICULTURAL LAND — DISTANCE TO FOREST RESERVE
################################################################################

## Standardise CRS: project forest reserves to UTM 36N
forest_reserves_utm <- forest_reserves |> st_transform(32636)

## ============================== Bushenyi ===================================

## Convert Bushenyi cropland points to sf and project
points_bushenyi_sf <- st_as_sf(points_bushenyi) |> st_transform(32636)

## Prepare observed distances for comparison
TIST_bushenyi <- TistDist_Dat |>
  dplyr::filter(Proj_Area == "Bushenyi") |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(32636)

nrow(TIST_bushenyi)

## Distance matrix (points x reserves) -> nearest distance (m)
dist_mat_rand <- st_distance(points_bushenyi_sf, forest_reserves_utm)
rand_dist_to_forest <- apply(dist_mat_rand, 1, min)
points_bushenyi_sf$Dist_To_Forest_m <- as.numeric(rand_dist_to_forest)

## Assemble comparison dataframe
obs_df <- TIST_bushenyi |>
  dplyr::select(Dist_To_Forest_m) |>
  dplyr::mutate(type = "TIST")

rand_df <- points_bushenyi_sf |>
  dplyr::select(Dist_To_Forest_m) |>
  dplyr::mutate(type = "Random_ag")

compare_df <- dplyr::bind_rows(obs_df, rand_df) |>
  dplyr::filter(!is.na(Dist_To_Forest_m),   # handles missing values
                Dist_To_Forest_m >= 0)

## ---- Diagnostic density plot ---------------------------------------------------

p_Bush_density <- ggplot(compare_df, aes(x = Dist_To_Forest_m, fill = type)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x     = "Distance to nearest forest reserve (m)",
    y     = "Density",
    fill  = NULL,
    title = "Bushenyi"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(size = 13, colour = "black"),
    axis.text.x      = element_text(angle = 90, hjust = 1, colour = "black"),
    axis.title       = element_text(size = 16, colour = "black"),
    legend.text      = element_text(size = 13, colour = "black"),
    legend.title     = element_text(size = 14, colour = "black")
  )
p_Bush_density

ggsave("Output/Manuscript 3 graphs/p_Bush_density.png", p_Bush_density, dpi = 600, width = 12, height = 6)

## ---- Compare distributions ------------------------------------------------------

compare_df |>
  sf::st_drop_geometry() |>
  dplyr::group_by(type) |>
  dplyr::summarise(
    q25 = quantile(Dist_To_Forest_m, 0.25, na.rm = TRUE),
    q50 = quantile(Dist_To_Forest_m, 0.50, na.rm = TRUE),
    q75 = quantile(Dist_To_Forest_m, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
## non-parametric statistical test
# type        q25   q50    q75
# <chr>     <dbl> <dbl>  <dbl>
# Random_ag 3647. 9424. 15408.
# TIST       831. 2874.  5945.

wilcox.test(
  x = compare_df$Dist_To_Forest_m[compare_df$type == "TIST"],
  y = compare_df$Dist_To_Forest_m[compare_df$type == "Random_ag"],
  alternative = "less"
)  # p-value < 0.05 (2.2e-16): TIST groves significantly closer to forest reserves than randomly sampled plots

## Compare median and mean to check credibility of Wilcoxon result
compare_df |>
  sf::st_drop_geometry() |>
  dplyr::group_by(type) |>
  dplyr::summarise(mean = mean(Dist_To_Forest_m), median = median(Dist_To_Forest_m))
# type       mean median
# <chr>     <dbl>  <dbl>
# Random_ag 9955.  9424.
# TIST      4477.  2874.

## ================================ Soroti ====================================

## Convert random cropland points to sf and project to UTM
points_soroti_sf <- st_as_sf(points_soroti) |> st_transform(32636)

## Prepare observed TIST groves for Soroti
TIST_soroti <- TistDist_Dat |>
  dplyr::filter(Proj_Area == "Soroti") |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(32636)

nrow(TIST_soroti)  # check sample size

## Distance matrix (points x forest reserves) -> nearest distance (m)
dist_mat_rand_soroti <- st_distance(points_soroti_sf, forest_reserves_utm)
rand_dist_to_forest_soroti <- apply(dist_mat_rand_soroti, 1, min)
points_soroti_sf$Dist_To_Forest_m <- as.numeric(rand_dist_to_forest_soroti)

## Assemble comparison dataframe
obs_df_soroti <- TIST_soroti |>
  dplyr::select(Dist_To_Forest_m) |>
  dplyr::mutate(type = "TIST")

rand_df_soroti <- points_soroti_sf |>
  dplyr::select(Dist_To_Forest_m) |>
  dplyr::mutate(type = "Random_ag")

compare_df_soroti <- dplyr::bind_rows(obs_df_soroti, rand_df_soroti) |>
  dplyr::filter(!is.na(Dist_To_Forest_m), Dist_To_Forest_m >= 0)

## ---- Diagnostic density plot ---------------------------------------------------

p_Soroti_density <- ggplot(compare_df_soroti, aes(x = Dist_To_Forest_m, fill = type)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x     = "Distance to nearest forest reserve (m)",
    y     = "Density",
    fill  = NULL,
    title = "Soroti"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(size = 13, colour = "black"),
    axis.text.x      = element_text(angle = 90, hjust = 1, colour = "black"),
    axis.title       = element_text(size = 16, colour = "black"),
    legend.text      = element_text(size = 13, colour = "black"),
    legend.title     = element_text(size = 14, colour = "black")
  )
p_Soroti_density

ggsave("Output/Manuscript 3 graphs/p_Soroti_density.png", p_Soroti_density, dpi = 600, width = 12, height = 6)


################################################################################
# 10. STANDARDISE Y-AXIS ACROSS BOTH DENSITY PLOTS
################################################################################
#
# geom_density() auto-scales y to each plot's own data, so comparing peak
# heights across Bushenyi vs Soroti panels is misleading otherwise -- a
# taller peak in one panel may just reflect a different y-axis range, not
# a real difference in concentration. Extracting the actual computed density
# values (via ggplot_build()) and applying one shared ylim to both plots
# fixes this without altering the underlying data.
################################################################################

bush_max_density   <- max(ggplot_build(p_Bush_density)$data[[1]]$y)
soroti_max_density <- max(ggplot_build(p_Soroti_density)$data[[1]]$y)

shared_ymax <- max(bush_max_density, soroti_max_density) * 1.05  # 5% headroom so the tallest peak doesn't touch the border

cat("Bushenyi max density:", round(bush_max_density, 6), "\n")
cat("Soroti max density:  ", round(soroti_max_density, 6), "\n")
cat("Shared y-axis max:   ", round(shared_ymax, 6), "\n")

## coord_cartesian() only zooms the view -- unlike scale_y_continuous(limits=...),
## it never re-filters the underlying data or distorts the density curve shape.
p_Bush_density <- p_Bush_density +
  coord_cartesian(ylim = c(0, shared_ymax))

p_Soroti_density <- p_Soroti_density +
  coord_cartesian(ylim = c(0, shared_ymax))

p_Bush_density
p_Soroti_density

## Re-save each individual plot with the standardised y-axis
ggsave("Output/Manuscript 3 graphs/p_Bush_density.png",  p_Bush_density,  dpi = 600, width = 12, height = 6)
ggsave("Output/Manuscript 3 graphs/p_Soroti_density.png", p_Soroti_density, dpi = 600, width = 12, height = 6)


################################################################################
# 11. COMBINED PANEL: BUSHENYI vs SOROTI, SIDE BY SIDE
################################################################################

p_density_combined <- (p_Bush_density | p_Soroti_density) +
  plot_layout(guides = "collect") +
  plot_annotation(
    # title = "TIST groves vs eligible planting sites: distance to nearest forest reserve"
  ) &
  theme(legend.position = "bottom")

p_density_combined

ggsave(
  "Output/Manuscript 3 graphs/p_density_combined.png",
  p_density_combined,
  dpi    = 600,
  width  = 16,
  height = 6
)


################################################################################
# 12. COMBINE RIDGELINE (A, TOP) + DENSITY COMPARISON (B, BOTTOM)
#     INTO ONE LABELLED PANEL
################################################################################

img_ridgeline <- image_read("Output/Manuscript 3 graphs/Ridge_dist.png")
img_density   <- image_read("Output/Manuscript 3 graphs/p_density_combined.png")

## Match widths before stacking, so neither panel looks stretched/squashed
## relative to the other -- image_append(stack=TRUE) requires equal widths.
target_width <- max(image_info(img_ridgeline)$width, image_info(img_density)$width)

img_ridgeline <- image_resize(img_ridgeline, paste0(target_width, "x"))
img_density   <- image_resize(img_density,   paste0(target_width, "x"))

## Panel labels -- same style/size as the earlier a/b/c/d variogram
## combination step, for visual consistency across the manuscript.
img_ridgeline <- image_annotate(img_ridgeline, "A", size = 140, location = "+20+20", weight = 700)
img_density   <- image_annotate(img_density,   "B", size = 140, location = "+20+20", weight = 700)

p_forest_proximity_combined <- image_append(
  c(img_ridgeline, img_density),
  stack = TRUE   # stack = TRUE means top-to-bottom; A (ridgeline) goes first
)

p_forest_proximity_combined

image_write(
  p_forest_proximity_combined,
  path = "Output/Manuscript 3 graphs/p_forest_proximity_combined.png"
)


################################################################################
# 13. SOROTI — DISTRIBUTION SUMMARY & SIGNIFICANCE TEST
################################################################################

## Summarise quantiles
compare_df_soroti |>
  sf::st_drop_geometry() |>
  dplyr::group_by(type) |>
  dplyr::summarise(
    q25 = quantile(Dist_To_Forest_m, 0.25, na.rm = TRUE),
    q50 = quantile(Dist_To_Forest_m, 0.50, na.rm = TRUE),
    q75 = quantile(Dist_To_Forest_m, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
# type        q25   q50    q75
# <chr>     <dbl> <dbl>  <dbl>
# Random_ag 3484. 6256. 10132.
# TIST      3783. 6772. 11085.

## Non-parametric test
wilcox.test(
  x = compare_df_soroti$Dist_To_Forest_m[compare_df_soroti$type == "TIST"],
  y = compare_df_soroti$Dist_To_Forest_m[compare_df_soroti$type == "Random_ag"],
  alternative = "less"
)  # p-value > 0.05 (0.9965)

## Check means/medians for sanity
compare_df_soroti |>
  sf::st_drop_geometry() |>
  dplyr::group_by(type) |>
  dplyr::summarise(mean = mean(Dist_To_Forest_m), median = median(Dist_To_Forest_m))
# type       mean median
# <chr>     <dbl>  <dbl>
# Random_ag 7407.  6256.
# TIST      7448.  6772.

################################################################################
# END OF SCRIPT
################################################################################
