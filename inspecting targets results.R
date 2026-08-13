library(targets)

# 1. Confirm ERGM model results completed without error
tar_meta(fields = c("name", "type", "error", "time")) |>
  dplyr::filter(grepl("^ergm|^site_models|^site_networks|^model_diagnostics", name))

# 2. Confirm overlap results completed without error (same check as before,
#    worth re-running since tar_outdated() flagged the general ERGM branch)
tar_meta(fields = c("name", "type", "error", "time")) |>
  dplyr::filter(grepl("overlap", name, ignore.case = TRUE))

# 3. Pull the actual values directly from the store -- this does NOT
#    rebuild anything, just reads what's already saved
ergm_results <- tar_read(ergm_results)
print(ergm_results, n = Inf)

final_ergm_table <- tar_read(final_ergm_table)
print(final_ergm_table, n = Inf)

overlap_results <- tar_read(overlap_results)
print(overlap_results, n = Inf)

final_overlap_table <- tar_read(final_overlap_table)
print(final_overlap_table, n = Inf)

################
# ## network results confirmed to be present. now running lmem 
# tar_make(names = c(
#   "Bushenyi_base_density_model", "Soroti_base_density_model",
#   "Bushenyi_H4_density_models", "Soroti_H4_density_models",
#   "Bushenyi_base_density_ML", "Soroti_base_density_ML",
#   "mixed_model_convergence"
# ))

# ============================================================
# MIXED-EFFECTS DENSITY MODELS -- confirm the RE-structure fix worked
# ============================================================

# The single most direct check: did "Warning" become "OK"?
mixed_model_convergence <- tar_read(mixed_model_convergence)
print(mixed_model_convergence, n = Inf)

# Full coefficient tables + convergence/singularity messages for the
# base density models -- summary() will show "singular fit" text
# directly if the problem is still there.
Bushenyi_base_density_model <- tar_read(Bushenyi_base_density_model)
summary(Bushenyi_base_density_model)

Soroti_base_density_model <- tar_read(Soroti_base_density_model)
summary(Soroti_base_density_model)

# ML refits (REML = FALSE) used for AIC comparison -- these inherit
# the fix automatically via update(), just confirming here.
Bushenyi_base_density_ML <- tar_read(Bushenyi_base_density_ML)
AIC(Bushenyi_base_density_ML)

Soroti_base_density_ML <- tar_read(Soroti_base_density_ML)
AIC(Soroti_base_density_ML)

# H4 density models -- each is a list of 10 fits (one per NearFar
# threshold). Pull AIC and the Exposure_sc coefficient across all
# 10 as a compact table, plus check the 1st threshold's full summary
# as a spot check for singularity.
Bushenyi_H4_density_models <- tar_read(Bushenyi_H4_density_models)
Soroti_H4_density_models   <- tar_read(Soroti_H4_density_models)

h4_density_check <- function(model_list, site_label) {
  tibble::tibble(
    site = site_label,
    threshold = seq_along(model_list),
    AIC = purrr::map_dbl(model_list, AIC),
    singular = purrr::map_lgl(model_list, lme4::isSingular, tol = 1e-4),
    Exposure_sc_estimate = purrr::map_dbl(
      model_list,
      ~ coef(summary(.x))["Exposure_sc", "Estimate"]
    )
  )
}

dplyr::bind_rows(
  h4_density_check(Bushenyi_H4_density_models, "Bushenyi"),
  h4_density_check(Soroti_H4_density_models, "Soroti")
) |> print(n = Inf)

# Spot check: full summary of Bushenyi's threshold-1 H4 density model
summary(Bushenyi_H4_density_models[[1]])

# ============================================================
# LAND-COVER
# ============================================================

# Per-site landcover class breakdown for TIST groves
Bushenyi_landcover <- tar_read(Bushenyi_landcover)
print(Bushenyi_landcover, n = Inf)

Soroti_landcover <- tar_read(Soroti_landcover)
print(Soroti_landcover, n = Inf)

# Suitability classification (eligible vs. excluded, by suitability class)
Bushenyi_suitability <- tar_read(Bushenyi_suitability)
print(Bushenyi_suitability, n = Inf)

Soroti_suitability <- tar_read(Soroti_suitability)
print(Soroti_suitability, n = Inf)

# Figures -- these are ggplot objects, printing renders them
print(tar_read(Bushenyi_landcover_plot))
print(tar_read(Soroti_landcover_plot))
print(tar_read(landcover_comparison_figure))
print(tar_read(Bushenyi_suitability_plot))
print(tar_read(Soroti_suitability_plot))


# ============================================================
# ELIGIBILITY / FOREST EXCLUSION
# These are format = "file" targets -- tar_read() gives you the
# path string; terra::rast() on top actually loads the raster.
# ============================================================

eligible_mask_file <- tar_read(eligible_mask_file)
cat("eligible_mask_file path:", eligible_mask_file, "\n")
eligible_mask <- terra::rast(eligible_mask_file)
terra::plot(eligible_mask, main = "Eligible landcover mask (pre-forest-exclusion)")

forest_mask_file <- tar_read(forest_mask_file)
forest_mask <- terra::rast(forest_mask_file)
terra::plot(forest_mask, main = "Rasterised forest reserves")

eligible_area_file <- tar_read(eligible_area_file)
eligible_area <- terra::rast(eligible_area_file)
terra::plot(eligible_area, main = "Eligible area (forest reserves excluded)")

eligible_Bushenyi_file <- tar_read(eligible_Bushenyi)
eligible_Bushenyi <- terra::rast(eligible_Bushenyi_file)
terra::plot(eligible_Bushenyi, main = "Eligible area -- Bushenyi")

eligible_Soroti_file <- tar_read(eligible_Soroti)
eligible_Soroti <- terra::rast(eligible_Soroti_file)
terra::plot(eligible_Soroti, main = "Eligible area -- Soroti")

# # Quick cell counts (eligible vs. NA), useful as a sanity check on
# # how much land actually survived the exclusion at each stage
# cat("Eligible cells (mask, pre-exclusion):", sum(!is.na(terra::values(eligible_mask))), "\n")
# cat("Eligible cells (post forest exclusion):", sum(!is.na(terra::values(eligible_area))), "\n")
# cat("Eligible cells -- Bushenyi:", sum(!is.na(terra::values(eligible_Bushenyi))), "\n")
# cat("Eligible cells -- Soroti:", sum(!is.na(terra::values(eligible_Soroti))), "\n")
# ============================================================
# FOREST-DISTANCE ANALYSIS
# ============================================================

# Summary statistics: TIST groves vs. random agricultural points
Bushenyi_distance_summary <- tar_read(Bushenyi_distance_summary)
print(Bushenyi_distance_summary, n = Inf)

Soroti_distance_summary <- tar_read(Soroti_distance_summary)
print(Soroti_distance_summary, n = Inf)

# Wilcoxon test results -- is TIST placement significantly closer
# to forest reserves than random agricultural land?
Bushenyi_distance_test <- tar_read(Bushenyi_distance_test)
print(Bushenyi_distance_test)

Soroti_distance_test <- tar_read(Soroti_distance_test)
print(Soroti_distance_test)

# Raw comparison data, in case you want a quick density plot --
# no dedicated plot target exists for this in the current pipeline,
# but the underlying data is right here
Bushenyi_distance_comparison <- tar_read(Bushenyi_distance_comparison)
Soroti_distance_comparison <- tar_read(Soroti_distance_comparison)

ggplot2::ggplot(Bushenyi_distance_comparison, ggplot2::aes(x = Dist_To_Forest_m, fill = type)) +
  ggplot2::geom_density(alpha = 0.4) +
  ggplot2::labs(title = "Bushenyi: distance to nearest forest reserve")

ggplot2::ggplot(Soroti_distance_comparison, ggplot2::aes(x = Dist_To_Forest_m, fill = type)) +
  ggplot2::geom_density(alpha = 0.4) +
  ggplot2::labs(title = "Soroti: distance to nearest forest reserve")
