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
## network results confirmed to be present.
tar_make(names = c(
  "Bushenyi_base_density_model", "Soroti_base_density_model",
  "Bushenyi_H4_density_models", "Soroti_H4_density_models",
  "Bushenyi_base_density_ML", "Soroti_base_density_ML",
  "mixed_model_convergence"
))

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
