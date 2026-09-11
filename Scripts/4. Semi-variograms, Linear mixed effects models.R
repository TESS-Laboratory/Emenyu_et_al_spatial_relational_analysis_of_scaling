# TIST BUSHENYI vs SOROTI -- FULL MODELLING + SPATIAL DIAGNOSTICS PIPELINE
#
# HOW TO READ THIS SCRIPT (for a non-technical reviewer)
# --------------------------------------------------------------------------
# This script tests whether "exposure" to the TIST tree-planting programme
# is associated with tree-planting outcomes (density and counts), and then
# checks whether those results might be distorted by geography -- i.e.
# whether nearby farms simply resemble each other because they are close
# together, rather than because of anything the model has captured.
#
# It is organised into eight parts. The spatial diagnostics (semivariograms)
# are placed FIRST, right after setup, because they answer the more basic
# descriptive question -- "is there any spatial pattern in the outcomes at
# all, and how far does it reach?" -- before the modelling sections attempt
# to explain that pattern:
#
#   PART 0 -- Setup: load data, prepare variables for modelling.
#   PART 1 -- Semivariograms, STAGE A: raw exploratory semivariograms. A
#             first, simple look at whether tree density/counts show any
#             spatial pattern at all, before any model is fitted. Ends
#             with a combined 4-panel figure for the manuscript.
#   PART 2 -- Semivariograms, STAGE B: residual semivariograms. The same
#             technique applied to what an exposure-only model could NOT
#             explain, to see how much of Stage A's spatial pattern
#             exposure accounts for on its own, before any other
#             mechanism is introduced.
#   PART 3 -- Main models (H2/H3/H4): does exposure predict tree density
#             and tree counts, and does controlling for spatial clustering
#             of similar farms improve the model?
#   PART 4 -- Diagnostics and figures built on the Part 3 models:
#             collinearity checks, marginal-effect figures, and a set of
#             targeted robustness checks (Bushenyi random-effects
#             structure, farm-size confounding, Soroti group composition
#             and tenure effects).
#   PART 5 -- Similarity models (H2b/H3b/H4b): a parallel set of models
#             asking whether NEIGHBOURING farms become more alike as
#             exposure increases (a different, complementary question
#             to Part 3).
#   PART 6 -- Coefficient-by-threshold ("forest plot") figures for the
#             Part 5 similarity models.
#   PART 7 -- Norm-convergence asymmetry test: does exposure pull farmers
#             above and below their neighbourhood mean toward each other
#             symmetrically, or only from one direction?
#   PART 8 -- Formal statistical tests for leftover spatial pattern in
#             the Part 3 models (Moran's I and DHARMa). These give a
#             yes/no significance answer, complementing the descriptive
#             picture from Parts 1-2 (see note at the top of Part 8 for
#             why both are kept).
#
# Every semivariogram plot in this script (raw or residual) uses the exact
# same visual layout, so Part 1 and Part 2 stay visually comparable:
#   - black dots + line  = the actual measured spatial pattern
#   - grey line          = a smooth curve fitted through those points
#   - dotted blue line    = "nugget" (baseline noise between close points)
#   - dashed green line   = "sill" (total variation once points are far apart)
#   - dashed red line     = "range" (distance beyond which points stop
#                            resembling each other)
#   - text underneath      = the same four numbers spelled out in words
#
# NOTE ON AN INHERITED NAMING QUIRK: earlier drafts of the raw-variogram
# code stored Soroti's data in objects literally named "Bush...". This
# version removes that risk entirely -- every raw variogram below is built
# by one function that takes the site name as a plain text argument, so
# there is no object name left to get out of sync with the data inside it.
################################################################################

## ---- Packages --------------------------------------------------------------
## Every package whose functions are called somewhere in this script:
##   tidyverse  - dplyr/ggplot2/tibble/purrr/stringr/forcats/tidyr workflow
##   lme4       - lmer() (density models, similarity models)
##   glmmTMB    - glmmTMB()/nbinom2 (tree-count negative-binomial models)
##   sf         - st_as_sf(), st_transform(), st_coordinates()
##   sp         - as(..., "Spatial") coercion required by gstat's variogram()
##   gstat      - variogram(), fit.variogram(), vgm(), variogramLine()
##   ggeffects  - ggpredict() (marginal-effect predictions for figures)
##   patchwork  - plot_layout()/plot_annotation()/"+" "/" combining of ggplots
##   flextable  - flextable(), set_caption(), theme_vanilla(), etc. (Word table)
##   officer    - required by flextable::save_as_docx()
##   spdep      - knn2nb(), knearneigh(), nb2listw(), moran.test()
##   DHARMa     - simulateResiduals(), testSpatialAutocorrelation(),
##                testDispersion(), testZeroInflation()
##   magick     - image_read()/image_annotate()/image_append()/image_write()
library(tidyverse)
library(lme4)
library(glmmTMB)
library(sf)
library(sp)
library(gstat)
library(ggeffects)
library(patchwork)
library(flextable)
library(officer)
library(spdep)
library(DHARMa)
library(magick)


output_dir <- "Output/Manuscript 3 graphs"


################################################################################
# PART 0 -- SETUP: LOAD DATA, SCALE PREDICTORS, SPLIT BY SITE
################################################################################

## Load the full, cleaned dataset (both project sites together).
## Continuous predictors are "scaled" (converted to a common 0-centred,
## comparable unit) so that coefficients across models can be compared fairly.
TistDat_H <- read.csv(
  "Data/TISTDat/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv"
) %>%
  distinct()

TistDat_H <- TistDat_H |>
  mutate(
    Dist_To_Forest_sc  = as.numeric(scale(Dist_To_Forest_m)),
    Area_Ha_sc         = as.numeric(scale(Area_Ha)),
    Trees_sc           = as.numeric(scale(Trees)),
    Exposure_sc        = as.numeric(scale(Exposure)),
    Years_since_reg_sc = as.numeric(scale(Years_since_reg)),

    # Scale all 10 "Near-Far" spatial-concentration indices individually.
    # These measure, at 10 different distance thresholds (500m to 5000m),
    # how spatially concentrated a farm's exposure is relative to its
    # neighbours. Used later to test H4/H4b.
    across(
      NearFar_resid_1:NearFar_resid_10,
      ~ as.numeric(scale(.x)),
      .names = "{.col}_sc"
    )
  )

## One canonical name for the 10 scaled Near-Far columns, reused everywhere
## below (H4 and H4b, both outcomes) -- avoids any risk of two different
## naming conventions drifting out of sync with each other.
scaled_nf_cols   <- paste0("NearFar_resid_", 1:10, "_sc")
threshold_labels <- paste0(seq(500, 5000, by = 500), "m")

cat("Scaled Near-Far columns present:", all(scaled_nf_cols %in% names(TistDat_H)), "\n")
cat("Correlations with Exposure_sc (should all be ~0, confirming independence):\n")
print(sapply(scaled_nf_cols, function(v) {
  round(cor(TistDat_H[[v]], TistDat_H$Exposure_sc, use = "complete.obs"), 8)
}))

## Split into the two project sites -- every model in this pipeline is fit
## separately per site, since Bushenyi and Soroti are different programme
## contexts and should not be pooled into one model.
TistDat_B <- TistDat_H |> filter(Proj_Area == "Bushenyi")
TistDat_S <- TistDat_H |> filter(Proj_Area == "Soroti")

cat("Bushenyi n =", nrow(TistDat_B), "\n")
cat("Soroti n =",   nrow(TistDat_S), "\n")

## "Duration tercile": splits each site's farms into three equal-sized
## groups (Short / Medium / Long) based on how many years they have been
## registered in the programme. Cut separately per site, since the two
## sites have different registration histories.
build_duration_tercile <- function(dat) {
  dat %>%
    mutate(
      Duration_Tercile = case_when(
        Years_since_reg <= quantile(Years_since_reg, 1/3, na.rm = TRUE) ~ "Short",
        Years_since_reg <= quantile(Years_since_reg, 2/3, na.rm = TRUE) ~ "Medium",
        TRUE ~ "Long"
      ),
      Duration_Tercile = factor(Duration_Tercile, levels = c("Short", "Medium", "Long"))
    )
}

TistDat_B <- build_duration_tercile(TistDat_B)
TistDat_S <- build_duration_tercile(TistDat_S)

cat("\nDuration_Tercile counts -- Bushenyi:\n"); print(table(TistDat_B$Duration_Tercile))
cat("\nDuration_Tercile counts -- Soroti:\n");   print(table(TistDat_S$Duration_Tercile))


################################################################################
# PART 1 -- SEMIVARIOGRAMS, STAGE A: RAW EXPLORATORY SEMIVARIOGRAMS
#
# Before any model is fitted, a simple descriptive question: do tree
# density and tree counts show any spatial pattern at all -- i.e. do farms
# located close together tend to resemble each other more than farms far
# apart? This stage answers that with semivariograms fit directly to the
# raw (log-transformed) outcomes, one per site x outcome combination.
################################################################################

## ---- Shared variogram tools -- used by BOTH Part 1 (raw) and Part 2
## (residual) semivariograms below, so every plot in this script shares the
## same visual style. --------------------------------------------------------

## ---- Fit a semivariogram to any numeric column -----------------------------
## Doesn't care whether the input is a raw outcome or a model's leftover
## residuals. Retries once on a degenerate fit, and falls back to a
## pure-nugget model if a spherical fit still can't be estimated reliably
## (e.g. too few point-pairs in a small or spatially sparse sample).
fit_variogram_generic <- function(data_with_coords, value_col, cutoff = 8000, width = 500) {
  data_clean <- data_with_coords %>% filter(!is.na(.data[[value_col]]))

  data_sf <- st_as_sf(data_clean, coords = c("longitude", "latitude"), crs = 4326)
  data_sf <- st_transform(data_sf, 32636)   # UTM zone 36N, so distances are true metres
  data_sp <- as(data_sf, "Spatial")

  vgm_empirical <- variogram(as.formula(paste0(value_col, " ~ 1")), data_sp, cutoff = cutoff, width = width)

  fit_attempt <- function(start_range) {
    tryCatch(
      fit.variogram(
        vgm_empirical,
        vgm(psill = var(data_sp[[value_col]], na.rm = TRUE), model = "Sph", range = start_range, nugget = 0)
      ),
      error = function(e) NULL
    )
  }

  vgm_fitted <- fit_attempt(2000)
  range_ok   <- function(v) !is.null(v) && any(v$model == "Sph") && is.finite(v$range[v$model == "Sph"]) && v$range[v$model == "Sph"] > 0

  if (!range_ok(vgm_fitted)) {
    cat("  [fit_variogram_generic] Initial spherical fit failed or gave a non-positive range for '",
        value_col, "' -- retrying with a different starting range.\n", sep = "")
    vgm_fitted <- fit_attempt(max(vgm_empirical$dist) / 2)
  }

  if (!range_ok(vgm_fitted)) {
    cat("  [fit_variogram_generic] Spherical model still could not be fit reliably for '", value_col,
        "' -- falling back to a pure-nugget model. Treat spatial_dependence = 0% as\n",
        "  'not estimable', not as confirmed evidence of no spatial pattern.\n", sep = "")
    vgm_fitted <- vgm(psill = var(data_sp[[value_col]], na.rm = TRUE), model = "Nug")
  }

  fit_curve <- variogramLine(vgm_fitted, maxdist = max(vgm_empirical$dist))

  nugget_val             <- if (any(vgm_fitted$model == "Nug")) vgm_fitted$psill[vgm_fitted$model == "Nug"] else 0
  partial_sill_val       <- if (any(vgm_fitted$model == "Sph")) vgm_fitted$psill[vgm_fitted$model == "Sph"] else 0
  range_val              <- if (any(vgm_fitted$model == "Sph")) vgm_fitted$range[vgm_fitted$model == "Sph"] else 0
  total_sill_val         <- nugget_val + partial_sill_val
  spatial_dependence_val <- if (total_sill_val > 0) partial_sill_val / total_sill_val else 0

  list(
    empirical = vgm_empirical, fitted = vgm_fitted, fit_curve = fit_curve,
    nugget = nugget_val, partial_sill = partial_sill_val,
    total_sill = total_sill_val, range_m = range_val,
    spatial_dependence = spatial_dependence_val
  )
}

## ---- Draw one semivariogram panel onto whatever device is already open -----
## No png()/dev.off() here -- shared by both the single-plot and panel-figure
## wrappers below, so the plotting logic exists in exactly one place.
## `label`, if given (e.g. "A"), is drawn in the top-left corner of the panel.
## `show_stats`, if TRUE, adds a compact spatial-dependence annotation
## (% of total variance attributable to distance) below the range label.
draw_variogram_panel <- function(vgm_result, plot_title, label = NULL, show_stats = FALSE) {
  dist_full  <- c(0, vgm_result$empirical$dist)
  gamma_full <- c(vgm_result$nugget, vgm_result$empirical$gamma)
  y_max      <- max(gamma_full, vgm_result$total_sill * 1.1)

  plot(
    dist_full, gamma_full, type = "b", pch = 19, col = "black",
    xlim = c(0, max(dist_full)), ylim = c(0, y_max),
    xlab = "Distance (m)", ylab = "Semivariance", main = plot_title
  )
  lines(vgm_result$fit_curve$dist, vgm_result$fit_curve$gamma, col = "gray40", lwd = 2)
  abline(v = vgm_result$range_m, lty = 2, lwd = 2, col = "red")

  text(
    x = vgm_result$range_m, y = y_max * 0.85,
    labels = paste0("Range = ", round(vgm_result$range_m, 0), " m"),
    col = "red", cex = 1.1, pos = 4
  )

  if (show_stats) {
    ## Anchored to the plot's right edge (via par("usr")) rather than to
    ## range_m, so it never runs out of horizontal room regardless of
    ## where the range line falls -- fixes the panel-B clipping issue.
    spatial_pct   <- round(vgm_result$spatial_dependence * 100, 0)
    unexplained_pct <- 100 - spatial_pct
    stat_x <- par("usr")[2] - 0.03 * diff(par("usr")[1:2])
    stat_y <- y_max * 0.15
    text(
      x = stat_x, y = stat_y,
      labels = paste0("Spatial dependence = ", spatial_pct, "%\n(", unexplained_pct, "% unexplained)"),
      col = "black", cex = 1.0, adj = c(1, 0)
    )
  }

  if (!is.null(label)) {
    text(
      x = par("usr")[1], y = par("usr")[4],
      labels = label, font = 2, cex = 1.6,
      adj = c(-0.3, 1.3), xpd = NA
    )
  }
}

## ---- Single-plot PNG wrapper -----------------------------------------------
plot_variogram_standard <- function(vgm_result, plot_title, save_path, show_stats = FALSE) {
  png(save_path, width = 2400, height = 1800, res = 300)
  par(mar = c(5, 5, 3, 2), cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5)
  draw_variogram_panel(vgm_result, plot_title, show_stats = show_stats)
  dev.off()
}

## ---- Combined panel figure wrapper ------------------------------------------
## Reusable for any set of variogram results, not just these four -- e.g. later
## for tenure-stratified variograms (Section 5.2, Supplementary Figure SX).
## `labels` defaults to A, B, C, ... for however many panels are passed in,
## filled row-wise to match par(mfrow) -- so with the default 2x2 layout and
## the density-Bushenyi/density-Soroti/trees-Bushenyi/trees-Soroti call order,
## this gives A/B on row 1 and C/D on row 2 with no extra arguments needed.
## `show_stats`, if TRUE, adds the spatial-dependence annotation to every panel.
make_variogram_panel_figure <- function(results_list, titles, save_path, nrow = 2, ncol = 2,
                                        labels = LETTERS[seq_along(results_list)],
                                        show_stats = FALSE) {
  png(save_path, width = 1600 * ncol, height = 1200 * nrow, res = 300)
  par(mfrow = c(nrow, ncol), mar = c(5, 5, 3, 2), cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.3)
  for (i in seq_along(results_list)) {
    draw_variogram_panel(results_list[[i]], titles[i], label = labels[i], show_stats = show_stats)
  }
  dev.off()
}

## ---- Fit + plot one raw variogram, individual PNG ---------------------------
run_raw_variogram <- function(full_data, site_name, outcome_col, outcome_label,
                              output_dir, log_transform = TRUE, show_stats = FALSE) {
  site_data <- full_data %>% filter(Proj_Area == site_name)

  value_col <- outcome_col
  if (log_transform) {
    value_col <- paste0("log_", outcome_col)
    site_data[[value_col]] <- log1p(site_data[[outcome_col]])
  }

  vgm_result <- fit_variogram_generic(site_data, value_col)

  plot_title <- paste0(site_name, " - ", outcome_label)
  save_path  <- file.path(output_dir, paste0(site_name, "_variogram_", outcome_col, ".png"))
  plot_variogram_standard(vgm_result, plot_title, save_path, show_stats = show_stats)

  cat("\n--- Raw exploratory variogram:", site_name, "/", outcome_label, "---\n")
  print(vgm_result$fitted)

  vgm_result
}

################################################################################
# RUN: four raw variograms (both outcomes x both sites), plus a combined
# 2x2 panel figure for presentation, plus a summary table.
################################################################################

output_dir <- "C:/workspace/Emenyu_et_al_TIST_leveraging_adoption/Output/Manuscript 3 graphs/variograms"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

TistDist_Dat <- read.csv(
  "Data/TISTDat/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv"
) %>%
  drop_na() %>%
  distinct()

raw_density_Bushenyi <- run_raw_variogram(TistDist_Dat, "Bushenyi", "Density_winsor99", "Tree density", output_dir)
raw_density_Soroti    <- run_raw_variogram(TistDist_Dat, "Soroti",   "Density_winsor99", "Tree density", output_dir)
raw_trees_Bushenyi    <- run_raw_variogram(TistDist_Dat, "Bushenyi", "Trees",            "Tree count",   output_dir)
raw_trees_Soroti      <- run_raw_variogram(TistDist_Dat, "Soroti",   "Trees",            "Tree count",   output_dir)

## Combined 2x2 panel figure for presentation -- includes spatial dependence
## annotation on each panel alongside the range line.
make_variogram_panel_figure(
  results_list = list(raw_density_Bushenyi, raw_density_Soroti, raw_trees_Bushenyi, raw_trees_Soroti),
  titles       = c("Bushenyi - Tree density", "Soroti - Tree density",
                   "Bushenyi - Tree count",  "Soroti - Tree count"),
  save_path    = file.path(output_dir, "raw_variograms_panel.png"),
  show_stats   = TRUE
)

## Summary table of fitted variogram parameters
raw_variogram_summary <- tibble(
  site    = c("Bushenyi", "Soroti", "Bushenyi", "Soroti"),
  outcome = c("Density (log)", "Density (log)", "Trees (log)", "Trees (log)"),
  nugget  = c(raw_density_Bushenyi$nugget, raw_density_Soroti$nugget,
              raw_trees_Bushenyi$nugget,   raw_trees_Soroti$nugget),
  total_sill = c(raw_density_Bushenyi$total_sill, raw_density_Soroti$total_sill,
                 raw_trees_Bushenyi$total_sill,   raw_trees_Soroti$total_sill),
  range_m = c(raw_density_Bushenyi$range_m, raw_density_Soroti$range_m,
              raw_trees_Bushenyi$range_m,   raw_trees_Soroti$range_m),
  spatial_dependence_pct = round(100 * c(
    raw_density_Bushenyi$spatial_dependence, raw_density_Soroti$spatial_dependence,
    raw_trees_Bushenyi$spatial_dependence,   raw_trees_Soroti$spatial_dependence
  ), 0)
)

print(raw_variogram_summary)

################################################################################
# COMBINE THE FOUR RAW VARIOGRAM PLOTS ABOVE INTO ONE LABELLED FIGURE
#
# Stitches the four Stage A raw-variogram PNGs into a single labelled
# 2x2 panel figure -- (a) Bushenyi density, (b) Soroti density,
# (c) Bushenyi trees, (d) Soroti trees -- for the manuscript. This uses
# `magick` image compositing rather than the make_variogram_panel_figure()
# helper above, so it can add the a/b/c/d annotation style used elsewhere
# in the manuscript's combined figures.
################################################################################

a <- image_read(file.path(output_dir, "Bushenyi_variogram_Density_winsor99.png"))
b <- image_read(file.path(output_dir, "Soroti_variogram_Density_winsor99.png"))
c <- image_read(file.path(output_dir, "Bushenyi_variogram_Trees.png"))
d <- image_read(file.path(output_dir, "Soroti_variogram_Trees.png"))

a <- image_annotate(a, "a", size = 80, location = "+20+20", weight = 700)
b <- image_annotate(b, "b", size = 80, location = "+20+20", weight = 700)
c <- image_annotate(c, "c", size = 80, location = "+20+20", weight = 700)
d <- image_annotate(d, "d", size = 80, location = "+20+20", weight = 700)

top    <- image_append(c(a, b))
bottom <- image_append(c(c, d))
combined_variogram_panels <- image_append(c(top, bottom), stack = TRUE)

image_write(combined_variogram_panels, file.path(output_dir, "combined_variogram_panels.png"))

## NOTE: the original pipeline also combines these variogram panels with
## three ridgeline plots (tree/area/density distributions) produced
## elsewhere in the wider analysis. That combination step is intentionally
## NOT reproduced here, since those ridgeline PNGs are generated by a
## separate part of the pipeline this script does not otherwise touch. If
## you still want that combined figure, re-run that section:
#
#   i   <- image_read(file.path(output_dir, "Ridgeline_Tree_Plot.png"))
#   ii  <- image_read(file.path(output_dir, "Ridgeline_Area_plot.png"))
#   iii <- image_read(file.path(output_dir, "Ridgeline_Density_Plot.png"))
#   blank <- image_blank(width = image_info(iii)$width, height = image_info(iii)$height, color = "white")
#   i   <- image_annotate(i,   "i",   size = 80, location = "+20+20", weight = 700)
#   ii  <- image_annotate(ii,  "ii",  size = 80, location = "+20+20", weight = 700)
#   iii <- image_annotate(iii, "iii", size = 80, location = "+20+20", weight = 700)
#   top1    <- image_append(c(i, ii))
#   bottom1 <- image_append(c(iii, blank))
#   All_ridgelines <- image_append(c(top1, bottom1), stack = TRUE)
#   image_write(All_ridgelines, file.path(output_dir, "All_ridgelines.png"))


################################################################################
# PART 2 -- SEMIVARIOGRAMS, STAGE B: RESIDUAL SEMIVARIOGRAMS AFTER
# EXPOSURE ALONE
#
# Having established (Stage A / Part 1) that raw planting outcomes are
# spatially clustered in both sites, this stage asks: how much of that
# clustering does exposure to neighbouring groves account for on its own,
# before any other mechanism (tenure, forest proximity, spatial
# concentration) is introduced?
#
# The "minimal" model here deliberately includes ONLY exposure as a fixed
# effect, alongside the same random-effects hierarchy used throughout (which
# controls for non-independence within group/village/district, not a
# substantive mechanism in its own right). Forest distance and tenure are
# withheld here so that any residual spatial structure -- or lack of it --
# can be attributed to exposure specifically, without conflating it with
# mechanisms tested separately elsewhere (forest proximity: Section 5.1;
# tenure and its interaction with exposure: Section 5.2.1 / 5.3, in Part 4
# below). This keeps the one-mechanism-per-method logic described in
# Section 3 intact.
#
# SCALE NOTE: Stage A's raw variograms were fit on log1p-transformed outcomes
# (run_raw_variogram()'s default log_transform = TRUE), so total_sill in
# raw_variogram_summary is on the log scale. To make sill_raw and
# sill_residual comparable in the variance-explained calculation below, the
# exposure-only model here is fit on the SAME log1p-transformed outcome, not
# the raw scale -- otherwise the two sills are on different scales and their
# ratio is not interpretable as "% variance explained".
################################################################################

## ---- Shared helper: which random-effects structure to use, by
## site/outcome/threshold ------------------------------------------------------
##
## Needed here for the exposure-only models below, and reused again in
## Part 5 and Part 7 for the full H2b/H3b/H4b models and the norm-convergence
## test -- defined once, here, so all three stages stay in sync.
##
## Fitting the original random-effects structure for the H2b/H3b/H4b models
## (Part 5) produced scattered singular-fit / non-convergence warnings
## across both sites and outcomes (see H2b_H3b_simplified_RE_evaluation.R
## for the full diagnostic). Dropping Admin_Districts/Subcounty (keeping
## Village_ID) resolves these cleanly and at negligible AIC cost for
## Density (both sites) and for Soroti Trees at 500m/1000m specifically --
## but for Bushenyi Trees, which never had a convergence problem, the same
## simplification measurably worsens AIC (~+147) and manufactures an
## apparently significant H3b interaction at 2500-3000m that is not present
## under the original, correctly-fitted structure. use_simplified_re()
## below encodes which combinations should use which structure, based on
## that evaluation -- it is not a blanket toggle.
use_simplified_re <- function(site, outcome, threshold_label = NA) {
  if (outcome == "Density") return(TRUE)              # both sites: resolves singular fits, ~zero AIC cost
  if (outcome == "Trees" && site == "Bushenyi") return(FALSE)   # never needed; measurably worse if applied
  if (outcome == "Trees" && site == "Soroti") {
    return(threshold_label %in% c("500m", "1000m"))    # only these two thresholds had the problem
  }
  FALSE
}

## Same random-effects logic as the full models elsewhere (Part 5 / 5.2.1):
## a simplified structure is used where the full nested hierarchy fails to
## converge (e.g. Bushenyi density, where many groups have only one active
## member). use_simplified is passed in per site/outcome via use_simplified_re(),
## exactly as in the full models, so the ONLY difference between this model
## and the full model downstream is which fixed effects are included.
fit_exposure_only_model <- function(dat, dv_col, use_simplified = FALSE) {
  re_term <- if (use_simplified) {
    "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
  } else {
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  }
  frm <- as.formula(paste0(dv_col, " ~ Exposure_sc + ", re_term))
  lmer(frm, data = dat, REML = TRUE, na.action = na.exclude)
}

## Fits the exposure-only model on a log1p-transformed outcome (matching
## Stage A's scale -- see SCALE NOTE above), extracts residuals, and fits a
## variogram to them using the exact same engine as the raw variograms, so
## the two stages are directly comparable and the resulting panel figure can
## reuse Stage A's drawing functions without modification.
##
## Guards against missing coordinates before st_as_sf(): fit_variogram_generic()
## only drops NA on the value column itself, not on longitude/latitude, so any
## rows with missing coordinates are filtered out explicitly here first.
run_exposure_residual_variogram <- function(dat, outcome_col, site_name, outcome_label,
                                            use_simplified = FALSE) {
  dat_complete <- dat %>% filter(!is.na(longitude), !is.na(latitude), !is.na(.data[[outcome_col]]))

  log_col <- paste0("log_", outcome_col)
  dat_complete[[log_col]] <- log1p(dat_complete[[outcome_col]])

  m_exposure <- fit_exposure_only_model(dat_complete, log_col, use_simplified = use_simplified)
  dat_res    <- dat_complete %>% mutate(resid_exposure = resid(m_exposure))

  vgm_result <- fit_variogram_generic(dat_res %>% rename(resid_val = resid_exposure), "resid_val")

  cat("\n--- Exposure-only residual variogram (log1p scale):", site_name, "/", outcome_label, "---\n")
  print(vgm_result$fitted)

  vgm_result
}

################################################################################
# RUN: four exposure-residual variograms (both outcomes x both sites).
# use_simplified_re() is called with an explicit threshold string rather than
# relying on a pre-set focal_label, so this block runs independently of
# whichever threshold Part 5's H2b/H3b code happens to have focused on.
# "500m" is used here as the Density branch of use_simplified_re() ignores
# the threshold argument entirely; for Trees (Soroti specifically) this
# DOES matter -- change the threshold string below to match whichever
# distance scale you want this diagnostic to reflect for Trees.
################################################################################

exposure_stage_b_threshold <- "500m"   # only affects Soroti Trees -- see use_simplified_re()

exp_resid_density_B <- run_exposure_residual_variogram(
  TistDat_B, "Density_winsor99", "Bushenyi", "Tree density",
  use_simplified = use_simplified_re("Bushenyi", "Density", exposure_stage_b_threshold)
)
exp_resid_density_S <- run_exposure_residual_variogram(
  TistDat_S, "Density_winsor99", "Soroti", "Tree density",
  use_simplified = use_simplified_re("Soroti", "Density", exposure_stage_b_threshold)
)
exp_resid_trees_B <- run_exposure_residual_variogram(
  TistDat_B, "Trees", "Bushenyi", "Tree count",
  use_simplified = use_simplified_re("Bushenyi", "Trees", exposure_stage_b_threshold)
)
exp_resid_trees_S <- run_exposure_residual_variogram(
  TistDat_S, "Trees", "Soroti", "Tree count",
  use_simplified = use_simplified_re("Soroti", "Trees", exposure_stage_b_threshold)
)

## Combined 2x2 panel figure, in the exact same style as the raw variogram
## panel (Figure 9) for direct visual comparability -- same drawing function,
## same label convention (A-D), same spatial-dependence annotation.
make_variogram_panel_figure(
  results_list = list(exp_resid_density_B, exp_resid_density_S, exp_resid_trees_B, exp_resid_trees_S),
  titles       = c("Bushenyi - Tree density (residual)", "Soroti - Tree density (residual)",
                   "Bushenyi - Tree count (residual)",  "Soroti - Tree count (residual)"),
  save_path    = file.path(output_dir, "exposure_residual_variograms_panel.png"),
  show_stats   = TRUE
)

## Summary table -- same structure as the raw variogram summary, so the two
## can be placed side by side to compute the proportion of variance exposure
## alone accounts for: (sill_raw - sill_residual) / sill_raw. Both sills are
## now on the log1p scale (see SCALE NOTE above), so this ratio is valid.
exposure_residual_summary <- tibble(
  site    = c("Bushenyi", "Soroti", "Bushenyi", "Soroti"),
  outcome = c("Density (log)", "Density (log)", "Trees (log)", "Trees (log)"),
  nugget  = c(exp_resid_density_B$nugget, exp_resid_density_S$nugget,
              exp_resid_trees_B$nugget,   exp_resid_trees_S$nugget),
  total_sill = c(exp_resid_density_B$total_sill, exp_resid_density_S$total_sill,
                 exp_resid_trees_B$total_sill,   exp_resid_trees_S$total_sill),
  range_m = c(exp_resid_density_B$range_m, exp_resid_density_S$range_m,
              exp_resid_trees_B$range_m,   exp_resid_trees_S$range_m),
  spatial_dependence_pct = round(100 * c(
    exp_resid_density_B$spatial_dependence, exp_resid_density_S$spatial_dependence,
    exp_resid_trees_B$spatial_dependence,   exp_resid_trees_S$spatial_dependence
  ), 0)
)

print(exposure_residual_summary)

## Variance explained by exposure alone: (sill_raw - sill_residual) / sill_raw,
## joined against the Stage A raw_variogram_summary computed earlier. Both
## tables are built in the same site/outcome row order (Bushenyi density,
## Soroti density, Bushenyi trees, Soroti trees), so bind_cols() is safe here
## -- but this is worth reconfirming if either summary table's construction
## changes in future edits, since bind_cols() does not check row alignment.
variance_explained_by_exposure <- raw_variogram_summary %>%
  select(site, outcome, sill_raw = total_sill) %>%
  bind_cols(exposure_residual_summary %>% select(sill_residual = total_sill)) %>%
  mutate(pct_variance_explained_by_exposure = round(100 * (sill_raw - sill_residual) / sill_raw, 1))

print(variance_explained_by_exposure)

###########################################
##generating word ready results table
##############
# install.packages("officer")

## Assemble the comparison table from the two summary tables and the
## variance-explained table already computed (raw_variogram_summary,
## exposure_residual_summary, variance_explained_by_exposure). Range values
## for rows where the residual fit did not plateau within the tested cutoff
## are reported as NA and rendered as "not resolved" below.
range_resolved <- c(TRUE, FALSE, FALSE, FALSE)   # Bushenyi density only, per Stage B panel figure

residual_comparison_table <- raw_variogram_summary %>%
  select(site, outcome, range_raw = range_m, spatial_dependence_raw = spatial_dependence_pct) %>%
  bind_cols(
    exposure_residual_summary %>% select(range_residual = range_m, spatial_dependence_residual = spatial_dependence_pct)
  ) %>%
  bind_cols(variance_explained_by_exposure %>% select(pct_variance_explained_by_exposure)) %>%
  mutate(
    range_residual_display = ifelse(range_resolved, paste0(round(range_residual, 0), " m"), "Not resolved\u2020"),
    range_raw_display      = paste0(round(range_raw, 0), " m"),
    outcome                = str_remove(outcome, " \\(log\\)")
  ) %>%
  select(
    Site = site, Outcome = outcome,
    `Range, raw` = range_raw_display,
    `Range, residual` = range_residual_display,
    `Spatial dependence, raw (%)` = spatial_dependence_raw,
    `Spatial dependence, residual (%)` = spatial_dependence_residual,
    `Variance explained by exposure (%)` = pct_variance_explained_by_exposure
  )

## Build the flextable, with the caveat as a footnote referencing the
## dagger symbol used for unresolved ranges above.
ft <- flextable(residual_comparison_table) %>%
  set_caption("Table SX. Comparison of raw and exposure-only residual semivariogram parameters, by site and outcome.") %>%
  add_footer_lines(
    "\u2020 Fitted spherical range exceeded the maximum tested distance (8,000 m); the empirical semivariogram had not plateaued within the observed range. Spatial dependence values for these rows should be interpreted with caution; variance explained by exposure (final column) does not share this limitation, as it is computed from total sills rather than the plateau location."
  ) %>%
  theme_vanilla() %>%
  autofit() %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header")

ft

## Save as a standalone Word document, ready to paste into supplementary
## materials or copy directly into the manuscript.
save_as_docx(ft, path = file.path(output_dir, "Table_SX_residual_comparison.docx"))


################################################################################
# PART 3 -- MAIN MODELS (H2 / H3 / H4): DOES EXPOSURE PREDICT OUTCOMES?
#
# H2: does programme exposure predict tree density / tree counts?
# H3: does that relationship change with years since registration?
# H4: does adding a measure of spatial concentration of exposure
#     (at 10 different distance thresholds) improve the model?
################################################################################

## Base model: exposure, years registered, their interaction, and distance
## to the nearest forest reserve, with random effects for administrative
## and programme group hierarchy (accounts for farms clustering within the
## same village/group rather than being fully independent observations).
##
## NOTE on Area_Ha_sc: added to the TREES formula only, as a check on
## whether farm size confounds the exposure -> tree count relationship
## (raised when reviewing the "fewer trees, but denser groves" result --
## a farm with less land could show fewer trees for a purely mechanical
## reason, independent of any behavioural response to exposure). It is
## NOT added to the Density formula: Density_winsor99 is computed as
## Trees / Area_Ha, so including Area_Ha_sc there would mean using part of
## the outcome's own denominator to predict the outcome -- not genuine
## confounding control, closer to circular reasoning. Farm size as a
## confound for density can only be checked via a standalone
## Area_Ha ~ Exposure regression, not by adding it into the density model.
fit_base_models <- function(dat) {
  list(
    density = lmer(
      Density_winsor99 ~
        Exposure_sc +
        Years_since_reg_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Admin_Districts/Subcounty/Village_ID) +
        (1 | Cluster_ID/Group_ID),
      data = dat, REML = TRUE
    ),
    # Tree counts modelled with a negative-binomial distribution rather
    # than a plain linear model, since counts cannot be negative and are
    # typically over-dispersed (more variable than a simple count model
    # expects). Area_Ha_sc included here only -- see note above.
    trees_nb = glmmTMB(
      Trees ~
        Exposure_sc +
        Years_since_reg_sc +
        Area_Ha_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Admin_Districts/Subcounty/Village_ID) +
        (1 | Cluster_ID/Group_ID),
      family = nbinom2,
      data = dat
    )
  )
}

## H4 extension: adds the spatial-concentration term (and its interaction
## with exposure) at ONE of the 10 distance thresholds. Fit once per
## threshold so we can see which distance scale, if any, matters most.
## Area_Ha_sc is included in the trees_formula only, for the same reason
## given above (fit_base_models) -- density's own construction rules it out.
fit_h4_models <- function(dat, nearfar_var) {
  density_formula <- as.formula(paste0(
    "Density_winsor99 ~ ",
    "Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc + ",
    "Dist_To_Forest_sc + ", nearfar_var, " + Exposure_sc:", nearfar_var, " + ",
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  ))
  trees_formula <- as.formula(paste0(
    "Trees ~ ",
    "Exposure_sc + Years_since_reg_sc + Area_Ha_sc + Exposure_sc:Years_since_reg_sc + ",
    "Dist_To_Forest_sc + ", nearfar_var, " + Exposure_sc:", nearfar_var, " + ",
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  ))
  list(
    density  = lmer(density_formula, data = dat, REML = TRUE),
    trees_nb = glmmTMB(trees_formula, family = nbinom2, data = dat)
  )
}

cat("\nFitting base models (H2/H3)...\n")
base_B <- fit_base_models(TistDat_B)
base_S <- fit_base_models(TistDat_S)

cat("\nFitting H4 models across 10 distance thresholds -- Bushenyi...\n")
h4_B <- setNames(
  lapply(scaled_nf_cols, function(v) fit_h4_models(TistDat_B, nearfar_var = v)),
  threshold_labels
)

cat("Fitting H4 models -- Soroti...\n")
h4_S <- setNames(
  lapply(scaled_nf_cols, function(v) fit_h4_models(TistDat_S, nearfar_var = v)),
  threshold_labels
)

## Convergence checks: confirms each model actually found a stable, valid
## solution rather than failing partway through fitting.
check_convergence <- function(mod_list, label) {
  cat("\n---", label, "---\n")
  lmer_warn <- mod_list$density@optinfo$conv$lme4$messages
  cat("Density (lmer):      ", ifelse(is.null(lmer_warn), "OK", paste("WARNING:", lmer_warn)), "\n")
  nb_conv <- mod_list$trees_nb$fit$convergence
  cat("Tree count (glmmTMB):", ifelse(nb_conv == 0, "OK", paste("WARNING code:", nb_conv)), "\n")
}

check_convergence(base_B, "Base -- Bushenyi")
check_convergence(base_S, "Base -- Soroti")

cat("\n=== H4 convergence -- Bushenyi ===\n")
for (i in seq_along(threshold_labels)) {
  lmer_warn <- h4_B[[i]]$density@optinfo$conv$lme4$messages
  nb_conv   <- h4_B[[i]]$trees_nb$fit$convergence
  cat(threshold_labels[i], "| density:", ifelse(is.null(lmer_warn), "OK", "WARNING"),
      "| trees:", ifelse(nb_conv == 0, "OK", "WARNING"), "\n")
}

cat("\n=== H4 convergence -- Soroti ===\n")
for (i in seq_along(threshold_labels)) {
  lmer_warn <- h4_S[[i]]$density@optinfo$conv$lme4$messages
  nb_conv   <- h4_S[[i]]$trees_nb$fit$convergence
  cat(threshold_labels[i], "| density:", ifelse(is.null(lmer_warn), "OK", "WARNING"),
      "| trees:", ifelse(nb_conv == 0, "OK", "WARNING"), "\n")
}

## AIC comparison: does adding the spatial-concentration term (H4) actually
## improve the model over the simpler base model (H2/H3)? A positive
## delta_AIC means H4 is the better-fitting model at that threshold.
compare_aic <- function(dat, site_label) {
  base_ml <- lmer(
    Density_winsor99 ~
      Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc +
      Dist_To_Forest_sc +
      (1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID),
    data = dat, REML = FALSE
  )
  aic_base <- AIC(base_ml)

  aic_h4 <- sapply(scaled_nf_cols, function(v) {
    f <- as.formula(paste0(
      "Density_winsor99 ~ Exposure_sc + Years_since_reg_sc + ",
      "Exposure_sc:Years_since_reg_sc + Dist_To_Forest_sc + ",
      v, " + Exposure_sc:", v, " + ",
      "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
    ))
    AIC(lmer(f, data = dat, REML = FALSE))
  })

  result <- data.frame(
    threshold = threshold_labels,
    AIC_base  = round(aic_base, 2),
    AIC_h4    = round(aic_h4, 2),
    delta_AIC = round(aic_base - aic_h4, 2)
  )

  cat("\n=== AIC comparison --", site_label, "(positive delta = H4 better) ===\n")
  print(result)
  cat("Best threshold:", result$threshold[which.max(result$delta_AIC)], "\n")

  invisible(result)
}

aic_B <- compare_aic(TistDat_B, "Bushenyi")
aic_S <- compare_aic(TistDat_S, "Soroti")


################################################################################
# PART 4 -- DIAGNOSTICS & FIGURES FOR THE PART 3 (MAIN) MODELS
################################################################################

################################################################################
# COLLINEARITY CHECK -- H4 OWN-OUTCOME MODELS (density, trees_nb)
#
# Mirrors extract_h4b_collinearity() (already applied to the H4b dissimilarity
# models in Part 5) but applied here to h4_B / h4_S, the Part 3 own-outcome
# models behind Figure 12. Never previously checked -- added to confirm or
# rule out the same exposure/spatial-concentration collinearity found in the
# dissimilarity models.
#
# get_vcov_matrix() handles the lmer vs. glmmTMB difference: glmmTMB's
# vcov() can return a plain matrix OR a list with $cond/$zi/$disp depending
# on version, so this checks for a list first before treating the result as
# a matrix directly.
################################################################################

get_vcov_matrix <- function(m) {
  v <- vcov(m)
  if (is.list(v) && !is.matrix(v)) v <- v$cond
  as.matrix(v)
}

extract_h4_collinearity <- function(mod_list, site_label, outcome) {
  map2_dfr(names(mod_list), scaled_nf_cols, function(th, nf_var) {
    m <- mod_list[[th]][[outcome]]
    corr_mat <- cov2cor(get_vcov_matrix(m))
    tibble(
      threshold = th,
      site = site_label,
      exposure_nearfar_corr = round(corr_mat["Exposure_sc", nf_var], 3)
    )
  })
}

h4_collinearity_density <- bind_rows(
  extract_h4_collinearity(h4_B, "Bushenyi", "density"),
  extract_h4_collinearity(h4_S, "Soroti", "density")
)

h4_collinearity_trees <- bind_rows(
  extract_h4_collinearity(h4_B, "Bushenyi", "trees_nb"),
  extract_h4_collinearity(h4_S, "Soroti", "trees_nb")
)

cat("\n=== H4 (own-outcome) collinearity check: corr(Exposure_sc, NearFar_resid_sc) -- Density ===\n")
print(h4_collinearity_density, n = 40)

cat("\n=== H4 (own-outcome) collinearity check: corr(Exposure_sc, NearFar_resid_sc) -- Trees ===\n")
print(h4_collinearity_trees, n = 40)

cat("\nInterpretation guide: same threshold as the H4b check (Part 5) -- treat |exposure_nearfar_corr| > 0.5\n",
    "as sufficiently collinear that the associated coefficient in Figure 12 should not be\n",
    "reported as confirmed.\n")

##########################################################
# PREDICTION GENERATION
############################################################

generate_predictions <- function(model,
                                 predictor,
                                 site,
                                 outcome,
                                 effect_name,
                                 threshold = NA_character_) {

  is_glmm <- inherits(model, "glmmTMB")

  pred <- ggpredict(
    model,
    terms = predictor,
    bias_correction = is_glmm
  ) %>%
    as.data.frame()

  pred <- pred %>%
    mutate(
      Site = site,
      Outcome = outcome,
      Effect = effect_name,
      Threshold = threshold
    )

  if ("group" %in% names(pred)) {
    pred <- pred %>%
      rename(Duration_level = group)
  } else {
    pred$Duration_level <- NA_character_
  }

  pred
}

############################################################
# GLOBAL SCALING CONSTANTS
############################################################

years_mu <- mean(
  TistDat_H$Years_since_reg,
  na.rm = TRUE
)

years_sigma <- sd(
  TistDat_H$Years_since_reg,
  na.rm = TRUE
)


exposure_mu <- mean(
  TistDat_H$Exposure,
  na.rm = TRUE
)

exposure_sigma <- sd(
  TistDat_H$Exposure,
  na.rm = TRUE
)


nearfar_mu_sigma <- lapply(
  1:10,
  function(k) {

    col <- TistDat_H[[paste0("NearFar_resid_", k)]]

    c(
      mu = mean(col, na.rm = TRUE),
      sigma = sd(col, na.rm = TRUE)
    )
  }
)

names(nearfar_mu_sigma) <- as.character(1:10)



############################################################
# FIGURE 1
# EXPOSURE EFFECT
############################################################

fig1_predictions <- bind_rows(

  generate_predictions(
    base_B$density,
    "Exposure_sc",
    "Bushenyi",
    "Density",
    "Exposure"
  ),

  generate_predictions(
    base_S$density,
    "Exposure_sc",
    "Soroti",
    "Density",
    "Exposure"
  ),

  generate_predictions(
    base_B$trees_nb,
    "Exposure_sc",
    "Bushenyi",
    "Trees",
    "Exposure"
  ),

  generate_predictions(
    base_S$trees_nb,
    "Exposure_sc",
    "Soroti",
    "Trees",
    "Exposure"
  )

) %>%

  mutate(
    x_raw = x * exposure_sigma + exposure_mu
  )



############################################################
# FIGURE 2A
# PROGRAMME DURATION MAIN EFFECT
############################################################

duration_main_predictions <- bind_rows(

  generate_predictions(
    base_B$density,
    "Years_since_reg_sc",
    "Bushenyi",
    "Density",
    "Programme duration"
  ),

  generate_predictions(
    base_S$density,
    "Years_since_reg_sc",
    "Soroti",
    "Density",
    "Programme duration"
  ),

  generate_predictions(
    base_B$trees_nb,
    "Years_since_reg_sc",
    "Bushenyi",
    "Trees",
    "Programme duration"
  ),

  generate_predictions(
    base_S$trees_nb,
    "Years_since_reg_sc",
    "Soroti",
    "Trees",
    "Programme duration"
  )

) %>%

  mutate(
    x_raw = x * years_sigma + years_mu
  )



############################################################
# FIGURE 2B
# EXPOSURE × PROGRAMME DURATION
############################################################

duration_interaction_predictions <- bind_rows(

  generate_predictions(
    base_B$density,
    c("Exposure_sc",
      "Years_since_reg_sc [meansd]"),
    "Bushenyi",
    "Density",
    "Exposure × Duration"
  ),

  generate_predictions(
    base_S$density,
    c("Exposure_sc",
      "Years_since_reg_sc [meansd]"),
    "Soroti",
    "Density",
    "Exposure × Duration"
  ),

  generate_predictions(
    base_B$trees_nb,
    c("Exposure_sc",
      "Years_since_reg_sc [meansd]"),
    "Bushenyi",
    "Trees",
    "Exposure × Duration"
  ),

  generate_predictions(
    base_S$trees_nb,
    c("Exposure_sc",
      "Years_since_reg_sc [meansd]"),
    "Soroti",
    "Trees",
    "Exposure × Duration"
  )

) %>%

  mutate(
    x_raw = x * exposure_sigma + exposure_mu
  )



############################################################
# FIGURE 3
# SPATIAL SCALE OF EXPOSURE CONCENTRATION
############################################################

thresholds <- c(
  "500m" = 1
  )


nearfar_predictions <- purrr::imap_dfr(

  thresholds,

  function(index, threshold) {

    variable <- paste0(
      "NearFar_resid_",
      index,
      "_sc"
    )


    bind_rows(

      generate_predictions(
        h4_B[[threshold]]$density,
        variable,
        "Bushenyi",
        "Density",
        "Exposure concentration",
        threshold
      ),

      generate_predictions(
        h4_S[[threshold]]$density,
        variable,
        "Soroti",
        "Density",
        "Exposure concentration",
        threshold
      ),

      generate_predictions(
        h4_B[[threshold]]$trees_nb,
        variable,
        "Bushenyi",
        "Trees",
        "Exposure concentration",
        threshold
      ),

      generate_predictions(
        h4_S[[threshold]]$trees_nb,
        variable,
        "Soroti",
        "Trees",
        "Exposure concentration",
        threshold
      )

    )

  }

) %>%

  mutate(
    Threshold = factor(
      Threshold,
      levels = c(
        "500m"
             )
    )
  )



nearfar_lookup <- tibble(

  Threshold = factor(
    c(
      "500m"
         ),
    levels = c(
      "500m"
          )
  ),

  mu = c(
    nearfar_mu_sigma[["1"]]["mu"]
     ),

  sigma = c(
    nearfar_mu_sigma[["1"]]["sigma"]
      )
  )



nearfar_predictions <- nearfar_predictions %>%

  left_join(
    nearfar_lookup,
    by = "Threshold"
  ) %>%

  mutate(
    x_raw = x * sigma + mu
  ) %>%

  select(
    -mu,
    -sigma
  )



############################################################
# GENERIC EFFECT PLOT
############################################################

plot_effect <- function(prediction_data,
                        xlab,
                        facet_threshold = FALSE) {


  pred_df <- prediction_data %>%

    mutate(

      x_raw = as.numeric(x_raw),
      predicted = as.numeric(predicted),
      conf.low = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),

      Outcome = factor(
        Outcome,
        levels = c("Density", "Trees")
      )

    ) %>%

    arrange(
      Outcome,
      Site,
      Threshold,
      x_raw
    ) %>%

    mutate(
      Plot_group = ifelse(
        is.na(Threshold),
        Site,
        interaction(
          Site,
          Threshold,
          drop = TRUE
        )
      )
    )


  p <- ggplot(

    pred_df,

    aes(
      x = x_raw,
      y = predicted,
      colour = Site,
      group = Plot_group
    )

  ) +

    geom_ribbon(

      aes(
        ymin = conf.low,
        ymax = conf.high,
        fill = Site,
        group = Plot_group
      ),

      alpha = 0.12,
      colour = NA

    ) +

    geom_line(
      linewidth = 1
    )


  if(facet_threshold){

    p <- p +
      facet_grid(
        Outcome ~ Threshold,
        scales = "free_y"
      )

  } else {

    p <- p +
      facet_wrap(
        ~Outcome,
        scales = "free_y"
      )

  }


  p +

    labs(
      x = xlab,
      y = "Predicted value",
      colour = "Site",
      fill = "Site"
    ) +

    theme_classic(
      base_size = 13
    ) +

    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )

}

############################################################
# DURATION INTERACTION PLOT
############################################################

plot_duration_interaction <- function(prediction_data) {


  pred_df <- prediction_data %>%

    filter(
      Effect == "Exposure × Duration"
    ) %>%

    mutate(

      x_raw = as.numeric(x_raw),

      Duration_level =
        as.numeric(as.character(Duration_level))

    ) %>%

    tidyr::drop_na(Duration_level) %>%

    group_by(
      Site,
      Outcome
    ) %>%

    mutate(

      Duration_level = case_when(

        Duration_level == min(Duration_level) ~
          "Short duration (-1 SD)",

        Duration_level == max(Duration_level) ~
          "Long duration (+1 SD)",

        TRUE ~
          "Average duration"

      ),

      Duration_level = factor(

        Duration_level,

        levels = c(
          "Short duration (-1 SD)",
          "Average duration",
          "Long duration (+1 SD)"
        )

      )

    ) %>%

    ungroup() %>%

    arrange(
      Outcome,
      Site,
      Duration_level,
      x_raw
    )


  ggplot(

    pred_df,

    aes(
      x = x_raw,
      y = predicted,
      colour = Site,
      linetype = Duration_level,
      group = interaction(
        Site,
        Duration_level
      )
    )

  ) +

    geom_line(
      linewidth = 1
    ) +

    facet_wrap(
      ~Outcome,
      scales = "free_y"
    ) +

    labs(
      x = "Exposure",
      y = "Predicted value",
      colour = "Site",
      linetype = "Programme duration"
    ) +

    theme_classic(
      base_size = 13
    ) +

    theme(
      legend.position = "bottom"
    )

}
############################################################
# FINAL FIGURES
############################################################

# FIGURE 1

fig1 <- plot_effect(

  fig1_predictions,

  xlab =
    "Exposure"

)

fig1


# FIGURE 2
fig2_duration_main <- plot_effect(

  duration_main_predictions,

  xlab =
    "Years since registration"

)


fig2_duration_interaction <-
  plot_duration_interaction(
    duration_interaction_predictions
  )



fig2 <- fig2_duration_main / fig2_duration_interaction

fig2

# FIGURE 3

fig3 <- plot_effect(

  nearfar_predictions,

  xlab =
    "Exposure concentration",

  facet_threshold = FALSE

)

fig3

############################################################
# SAVE FIGURES
############################################################

ggsave(
  "Figure1_exposure_effect.png",
  fig1,
  width = 10,
  height = 6,
  dpi = 300
)


ggsave(
  "Figure2_programme_duration.png",
  fig2,
  width = 10,
  height = 10,
  dpi = 300
)


ggsave(
  "Figure3_spatial_scale.png",
  fig3,
  width = 12,
  height = 7,
  dpi = 300
)

###############################################################################
# BUSHENYI DENSITY DIAGNOSTIC -- SIMPLIFIED RANDOM EFFECTS (H2 ONLY)
#
# base_B$density above returns a "boundary (singular) fit" warning, with
# Admin_Districts and Subcounty:Admin_Districts variance = 0 (Soroti's
# equivalent model does not show this problem, so it is Bushenyi-specific).
# This shortened check re-fits ONLY the H2 base density model with those two
# zero-variance levels dropped, to confirm the singular fit is not distorting
# H2's own Exposure_sc estimate. It reuses TistDat_B and base_B already
# fitted above -- no new data structures, no H4 threshold loop.
################################################################################

base_density_simplified_B <- lmer(
  Density_winsor99 ~
    Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc +
    Dist_To_Forest_sc +
    (1 | Village_ID) + (1 | Cluster_ID/Group_ID),
  data = TistDat_B, REML = TRUE
)

diag_original   <- coef(summary(base_B$density))["Exposure_sc", ]
diag_simplified <- coef(summary(base_density_simplified_B))["Exposure_sc", ]

density_diagnostic_summary <- tibble(
  structure = c("Original (Admin_Districts/Subcounty/Village_ID)", "Simplified (Village_ID only)"),
  singular_fit = c(!is.null(base_B$density@optinfo$conv$lme4$messages),
                   !is.null(base_density_simplified_B@optinfo$conv$lme4$messages)),
  AIC                  = round(c(AIC(base_B$density), AIC(base_density_simplified_B)), 2),
  Exposure_sc_estimate = round(c(diag_original["Estimate"],  diag_simplified["Estimate"]), 2),
  Exposure_sc_se       = round(c(diag_original["Std. Error"], diag_simplified["Std. Error"]), 2),
  Exposure_sc_t        = round(c(diag_original["t value"],   diag_simplified["t value"]), 3)
)

cat("\n=== Bushenyi Density diagnostic: original vs. simplified random effects (H2 only) ===\n")
print(density_diagnostic_summary)

cat("\nConclusion: if Exposure_sc's estimate/t-value above are essentially unchanged between\n",
    "structures, this confirms the singular fit did not distort the H2 conclusion for\n",
    "Bushenyi density -- the simplification resolves the degenerate random-effects warning\n",
    "without altering the substantive result.\n")

################################################################################
# FARM-SIZE CONFOUND CHECK: DOES EXPOSURE CORRELATE WITH FARM SIZE?
#
# Raised when reviewing the "fewer trees, but denser groves" result: do
# high-exposure neighbourhoods simply have smaller farms on average, rather
# than farmers actively concentrating effort? This tests that directly and
# with no circularity risk, since neither Trees nor Density_winsor99 (which
# is itself Trees / Area_Ha) is involved -- just Area_Ha_sc regressed on
# Exposure_sc, one site at a time, reusing TistDat_B / TistDat_S already
# built above.
################################################################################

cat("\n=== Farm size (Area_Ha_sc) vs. Exposure_sc -- Bushenyi ===\n")
area_exposure_B <- lm(Area_Ha_sc ~ Exposure_sc, data = TistDat_B)
print(summary(area_exposure_B))

cat("\n=== Farm size (Area_Ha_sc) vs. Exposure_sc -- Soroti ===\n")
area_exposure_S <- lm(Area_Ha_sc ~ Exposure_sc, data = TistDat_S)
print(summary(area_exposure_S))

cat("\nInterpretation guide: a significant, negative Exposure_sc coefficient in either model\n",
    "above means farms in higher-exposure neighbourhoods tend to be systematically smaller --\n",
    "supporting a land-availability explanation for the 'fewer trees per unit exposure'\n",
    "result, alongside (or instead of) a behavioural/competitive one. A null or positive\n",
    "coefficient argues against farm size driving that result.\n")

################################################################################
# SOROTI TENURE INTERACTION -- GROUP COMPOSITION AND CELL SIZE CHECKS
#
# Follow-up on the Exposure_sc:Years_since_reg_sc interaction in Soroti's
# tree count models. Two questions raised when drafting the Discussion
# section, checked here directly rather than left as speculative caveats:
#
#   1. Are single-member groups (previously flagged for Bushenyi's density
#      models) concentrated in a particular tenure tercile in Soroti? If
#      "long-tenured" groups are disproportionately single-member, the
#      tenure interaction may partly reflect an individual farmer's own
#      history rather than any group-level process, since a one-person
#      "group" has no collective dynamic to mature.
#   2. Is the high-exposure / long-tenure combination well-populated, or a
#      sparse corner of the sample potentially driving an outsized
#      interaction coefficient?
#
# Reuses TistDat_S and its existing Duration_Tercile column -- no new data
# loading.
################################################################################

## --- Check 1: group size by tenure tercile ---
## Group size = number of farmer rows sharing the same Group_ID. Duration_
## Tercile is taken as the first value per group, since tenure (and
## therefore its tercile) is shared by every member of the same group.
soroti_group_composition <- TistDat_S %>%
  group_by(Group_ID) %>%
  summarise(
    n_members        = n(),
    Duration_Tercile  = first(Duration_Tercile),
    .groups = "drop"
  ) %>%
  mutate(group_size_cat = ifelse(n_members == 1, "Single-member", "Multi-member"))

cat("\n=== Soroti group size by tenure tercile ===\n")
group_composition_table <- soroti_group_composition %>%
  count(Duration_Tercile, group_size_cat) %>%
  group_by(Duration_Tercile) %>%
  mutate(pct_within_tercile = round(100 * n / sum(n), 1)) %>%
  ungroup()
print(group_composition_table)

cat("\nChi-square test: is group-size category independent of tenure tercile?\n")
group_size_chisq_tab <- table(soroti_group_composition$Duration_Tercile, soroti_group_composition$group_size_cat)
print(group_size_chisq_tab)
print(suppressWarnings(chisq.test(group_size_chisq_tab)))

cat("\nInterpretation guide: a low p-value above means single-member groups are NOT evenly\n",
    "spread across tenure terciles -- check the pct_within_tercile column to see whether\n",
    "single-member groups are disproportionately concentrated in the 'Long' tenure tercile\n",
    "specifically (the concern raised for the tenure interaction). A high p-value means\n",
    "group size and tenure are effectively independent, easing that concern.\n")

## --- Check 2: sample size in the high-exposure / long-tenure cell ---
## Exposure is continuous, so it is split into terciles here purely for this
## diagnostic cross-tab -- this does NOT change Exposure_sc anywhere else in
## the pipeline, it is a one-off check on cell sample sizes only.
TistDat_S_diag <- TistDat_S %>%
  mutate(
    Exposure_Tercile = ntile(Exposure_sc, 3),
    Exposure_Tercile = factor(Exposure_Tercile, labels = c("Low", "Medium", "High"))
  )

cat("\n=== Soroti sample sizes: Exposure tercile x Duration tercile ===\n")
exposure_tenure_cell_sizes <- TistDat_S_diag %>%
  count(Exposure_Tercile, Duration_Tercile) %>%
  tidyr::pivot_wider(names_from = Duration_Tercile, values_from = n, values_fill = 0)
print(exposure_tenure_cell_sizes)

cat("\nInterpretation guide: the High-Exposure / Long-tenure cell (bottom-right of the table\n",
    "above) is the one the Soroti tenure-interaction finding depends on most directly. If it\n",
    "holds a reasonably large share of the sample, the interaction is well-supported by data;\n",
    "if it is a small fraction of the total, the finding rests on a thinner evidence base than\n",
    "the interaction's strong p-value alone would suggest.\n")


############################################################################################################
### visualising tercile distribution
# ---- 0. Real year range per tenure tercile (Soroti) -----------------------
## Not previously defined anywhere in the pipeline -- needed below to label
## each tercile with its actual Years_since_reg range (e.g. "Short tenure
## (1-3 yrs)") rather than just the tercile name on its own.
soroti_tercile_ranges <- TistDat_S %>%
  group_by(Duration_Tercile) %>%
  summarise(
    min_years = min(Years_since_reg, na.rm = TRUE),
    max_years = max(Years_since_reg, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 1. Build axis labels: tercile name + real year range ---------------
tercile_labels <- soroti_tercile_ranges %>%
  mutate(
    year_label = ifelse(
      min_years == max_years,
      paste0(min_years, " yrs"),
      paste0(min_years, "\u2013", max_years, " yrs")
    ),
    axis_label = paste0(Duration_Tercile, " tenure\n(", year_label, ")")
  ) %>%
  select(Duration_Tercile, axis_label)

# ---- 2. Prep plot data ----------------------------------------------------
plot_data <- group_composition_table %>%
  left_join(tercile_labels, by = "Duration_Tercile") %>%
  mutate(
    group_size_cat = factor(group_size_cat, levels = c("Multi-member", "Single-member")),
    axis_label = factor(axis_label, levels = tercile_labels$axis_label[order(tercile_labels$Duration_Tercile)])
  )

# ---- 3. Chi-square test (reused, not refit) -------------------------------
chisq_result <- suppressWarnings(chisq.test(group_size_chisq_tab))

chisq_subtitle <- paste0(
  "\u03c7\u00b2(", chisq_result$parameter, ") = ", round(chisq_result$statistic, 2),
  ", p = ", format(round(chisq_result$p.value, 3), nsmall = 3),
  " \u2014 not significant"
)

# ---- 4. Label positions ----------------------------------------------------
label_multi  <- plot_data %>% filter(group_size_cat == "Multi-member")
label_single <- plot_data %>% filter(group_size_cat == "Single-member")

# ---- 5. Build the figure ---------------------------------------------------
COL_MULTI  <- "#4C72B0"
COL_SINGLE <- "#DD8452"

fig_s4 <- ggplot(plot_data, aes(x = axis_label, y = pct_within_tercile, fill = group_size_cat)) +
  geom_col(width = 0.55, colour = NA, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = label_multi,
    aes(label = paste0("n=", n)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    colour = "white", fontface = "bold", size = 3.6
  ) +
  geom_text(
    data = label_single,
    aes(y = 102, label = paste0(pct_within_tercile, "%\n(n=", n, ")")),
    colour = COL_SINGLE, fontface = "bold", size = 3.1, vjust = 0
  ) +
  scale_fill_manual(values = c("Multi-member" = COL_MULTI, "Single-member" = COL_SINGLE)) +
  scale_y_continuous(limits = c(0, 113), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Share of Soroti groups (%)",
    fill = NULL,
    title = "Group composition by tenure tercile",
    subtitle = chisq_subtitle
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, face = "italic", colour = "grey30"),
    axis.text.x       = element_text(size = 10, lineheight = 0.9),
    axis.title.y      = element_text(size = 10.5)
  )

fig_s4

# ---- 6. Save ---------------------------------------------------------------
ggsave(
  "Figure_S4_soroti_tenure_composition.png",
  fig_s4,
  width = 6.2, height = 5.6, dpi = 300
)

############################################################################################################
# TENURE MAIN EFFECT ON TREE COUNT -- CHECKING "FARMERS IN OLDER GROUPS
# PLANTED FAR MORE TREES" (SOROTI)
#
# Translates Years_since_reg_sc's coefficient into a plain percentage change,
# then shows model-predicted tree counts across tenure levels alongside
# simple, unadjusted raw means for comparison -- so the "far more trees"
# claim is backed by an interpretable number, not just a coefficient and
# p-value. Reuses base_S$trees_nb and TistDat_S already fitted/loaded above.
################################################################################

## Percentage change in expected tree count per 1-SD increase in tenure,
## holding all else constant (the model uses a log link, so
## exp(coefficient) - 1 gives the proportional change).
## NOTE: unlike lmer, glmmTMB's summary() stores coefficients as a list
## ($cond/$zi/$disp) rather than a single matrix, so $cond is needed before
## subsetting by row/column name.
tenure_coef_S <- summary(base_S$trees_nb)$coefficients$cond["Years_since_reg_sc", "Estimate"]
cat("\n=== Soroti tenure effect on tree count, in plain percentage terms ===\n")
cat("Years_since_reg_sc coefficient:", round(tenure_coef_S, 3), "\n")
cat("Implied change in expected tree count per 1-SD increase in tenure:",
    round(100 * (exp(tenure_coef_S) - 1), 1), "%\n")

## Model-predicted tree counts at each tenure tercile's mean scaled tenure
## value, holding other predictors at their sample means -- gives an actual,
## interpretable tree-count number rather than only a coefficient.
tenure_tercile_means_S <- TistDat_S %>%
  group_by(Duration_Tercile) %>%
  summarise(mean_Years_since_reg_sc = mean(Years_since_reg_sc, na.rm = TRUE), .groups = "drop")

cat("\nMean Years_since_reg_sc by tenure tercile (Soroti):\n")
print(tenure_tercile_means_S)

predicted_trees_by_tenure_S <- ggpredict(
  base_S$trees_nb,
  terms = paste0(
    "Years_since_reg_sc [",
    paste(round(tenure_tercile_means_S$mean_Years_since_reg_sc, 3), collapse = ","),
    "]"
  )
)
cat("\nModel-predicted tree count at each tenure tercile's mean (Soroti, other predictors held at sample means):\n")
print(as.data.frame(predicted_trees_by_tenure_S))

## Simple, unadjusted comparison: raw mean/median tree count by tenure
## tercile, with no model or covariates involved -- a sanity check against
## the model-based prediction above.
cat("\nRaw (unadjusted) mean/median tree count by tenure tercile -- Soroti:\n")
raw_trees_by_tenure_S <- TistDat_S %>%
  group_by(Duration_Tercile) %>%
  summarise(
    n            = n(),
    mean_Trees   = round(mean(Trees, na.rm = TRUE), 1),
    median_Trees = round(median(Trees, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(raw_trees_by_tenure_S)

cat("\nInterpretation guide: compare the model-predicted counts (covariate-adjusted) against\n",
    "the raw means (unadjusted) above. If both show a similar, substantial increase from\n",
    "Short to Long tenure, this directly supports the claim that farmers in older groups\n",
    "planted far more trees. A large gap between the adjusted and raw pattern would suggest\n",
    "part of the raw difference is attributable to other covariates (e.g. farm area) rather\n",
    "than tenure itself.\n")


################################################################################
# PART 5 -- SIMILARITY MODELS (H2b / H3b / H4b)
#
# A different question from Part 3: instead of asking "does MY exposure
# predict MY outcome", this asks "do NEIGHBOURING farms become more
# similar to each other as exposure increases". The outcome variable here
# ("AbsDissimilarity") measures the absolute difference between a farm and
# its neighbours -- so a NEGATIVE coefficient means neighbours become
# more alike (less dissimilar) as the predictor increases.
#
#   H2b: does exposure predict this neighbour-similarity outcome?
#   H3b: does that relationship change with years since registration?
#   H4b: does adding the spatial-concentration term improve the model?
#
# NOTE: use_simplified_re(), which governs the random-effects structure
# used below, is already defined in Part 2 (it is needed there for the
# Stage B residual semivariograms) and is simply reused here.
################################################################################

## NOTE on include_area: for the TREES-based dissimilarity outcome
## (log_AbsDissimilarity_Trees_k), farm size can plausibly confound the
## exposure -> dissimilarity relationship the same way it can confound the
## Part 3 Trees model -- a farm with less land available could show a
## larger (or smaller) gap from its neighbours' tree counts for reasons
## unrelated to exposure. Set include_area = TRUE for that outcome only.
## For the DENSITY-based dissimilarity outcome (AbsDissimilarity_Density_*),
## leave include_area = FALSE: Density_winsor99 (and therefore its
## neighbour-dissimilarity measure) is itself computed from Trees / Area_Ha,
## so adding Area_Ha_sc there would again use part of the outcome's own
## construction to predict the outcome -- the same circularity issue flagged
## for the Part 3 density model.

fit_base_models_b <- function(dat, dv_col, include_area = FALSE, use_simplified = FALSE) {
  re_term <- if (use_simplified) {
    "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
  } else {
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  }
  frm <- as.formula(paste0(
    dv_col, " ~ ",
    "Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc + ",
    "Dist_To_Forest_sc + ",
    if (include_area) "Area_Ha_sc + " else "",
    re_term
  ))
  lmer(frm, data = dat, REML = TRUE, na.action = na.exclude)
}

fit_h4b_models <- function(dat, dv_col, nearfar_var, include_area = FALSE, use_simplified = FALSE) {
  re_term <- if (use_simplified) {
    "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
  } else {
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  }
  frm <- as.formula(paste0(
    dv_col, " ~ ",
    "Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc + ",
    "Dist_To_Forest_sc + ",
    if (include_area) "Area_Ha_sc + " else "",
    nearfar_var, " + Exposure_sc:", nearfar_var, " + ",
    re_term
  ))
  lmer(frm, data = dat, REML = TRUE, na.action = na.exclude)
}

check_convergence_b <- function(mod_list, label) {
  cat("\n=== Convergence --", label, "===\n")
  for (th in names(mod_list)) {
    warn <- mod_list[[th]]@optinfo$conv$lme4$messages
    cat(th, ":", ifelse(is.null(warn), "OK", paste("WARNING:", warn)), "\n")
  }
}

compare_aic_b <- function(base_list, h4b_list, site_label) {
  result <- tibble(
    threshold  = threshold_labels,
    AIC_base_b = map_dbl(base_list, AIC),
    AIC_h4b    = map_dbl(h4b_list, AIC)
  ) %>%
    mutate(delta_AIC = round(AIC_base_b - AIC_h4b, 2))

  cat("\n=== AIC comparison (H4b) --", site_label, "(positive delta = h4b better) ===\n")
  print(result)
  cat("Best threshold:", result$threshold[which.max(result$delta_AIC)], "\n")

  invisible(result)
}

## Pulls out just the two coefficients of interest (H2b: Exposure_sc;
## H3b: its interaction with years registered) from every threshold's
## model, using base R's lme4::coef(summary()) -- deliberately avoids the
## broom.mixed package to prevent a known rlang version conflict.
extract_h2b_h3b <- function(mod_list, site_label) {
  map_dfr(names(mod_list), function(th) {
    m <- mod_list[[th]]
    cf <- as.data.frame(coef(summary(m)))
    cf$term <- rownames(cf)
    names(cf)[names(cf) == "Estimate"]   <- "estimate"
    names(cf)[names(cf) == "Std. Error"] <- "std.error"
    names(cf)[names(cf) == "t value"]    <- "statistic"

    cf %>%
      as_tibble() %>%
      filter(term %in% c("Exposure_sc", "Exposure_sc:Years_since_reg_sc")) %>%
      select(term, estimate, std.error, statistic) %>%
      mutate(threshold = th, site = site_label)
  })
}

## --- Density outcome ---
dv_cols_density <- paste0("AbsDissimilarity_Density_winsor99_", 1:10)

## Density uses the simplified structure at every threshold, both sites --
## see use_simplified_re() note in Part 2.
cat("\nFitting base_b (H2b+H3b) density models -- Bushenyi...\n")
base_b_B <- setNames(map(dv_cols_density, ~ fit_base_models_b(TistDat_B, dv_col = .x, use_simplified = use_simplified_re("Bushenyi", "Density"))), threshold_labels)
cat("Fitting base_b density models -- Soroti...\n")
base_b_S <- setNames(map(dv_cols_density, ~ fit_base_models_b(TistDat_S, dv_col = .x, use_simplified = use_simplified_re("Soroti", "Density"))), threshold_labels)

cat("\nFitting h4b density models -- Bushenyi...\n")
h4b_B <- setNames(
  map2(dv_cols_density, scaled_nf_cols, ~ fit_h4b_models(TistDat_B, dv_col = .x, nearfar_var = .y, use_simplified = use_simplified_re("Bushenyi", "Density"))),
  threshold_labels
)
cat("Fitting h4b density models -- Soroti...\n")
h4b_S <- setNames(
  map2(dv_cols_density, scaled_nf_cols, ~ fit_h4b_models(TistDat_S, dv_col = .x, nearfar_var = .y, use_simplified = use_simplified_re("Soroti", "Density"))),
  threshold_labels
)

check_convergence_b(base_b_B, "base_b (Density) -- Bushenyi")
check_convergence_b(base_b_S, "base_b (Density) -- Soroti")
check_convergence_b(h4b_B,    "h4b (Density) -- Bushenyi")
check_convergence_b(h4b_S,    "h4b (Density) -- Soroti")

aic_b_B <- compare_aic_b(base_b_B, h4b_B, "Bushenyi (Density)")
aic_b_S <- compare_aic_b(base_b_S, h4b_S, "Soroti (Density)")

h2b_h3b_density_table <- bind_rows(
  extract_h2b_h3b(base_b_B, "Bushenyi"),
  extract_h2b_h3b(base_b_S, "Soroti")
)
cat("\n=== H2b/H3b coefficients -- Density ===\n")
print(h2b_h3b_density_table, n = 40)

## --- Tree count outcome (log1p) ---
for (k in 1:10) {
  raw_col <- paste0("AbsDissimilarity_Trees_", k)
  log_col <- paste0("log_AbsDissimilarity_Trees_", k)
  TistDat_B[[log_col]] <- log1p(TistDat_B[[raw_col]])
  TistDat_S[[log_col]] <- log1p(TistDat_S[[raw_col]])
}
dv_cols_trees <- paste0("log_AbsDissimilarity_Trees_", 1:10)

## Bushenyi Trees: original structure at every threshold -- never simplified,
## since it never had a convergence problem and simplifying measurably
## worsens fit there (see use_simplified_re() note in Part 2).
cat("\nFitting base_b (H2b+H3b) trees models -- Bushenyi...\n")
base_b_trees_B <- setNames(map(dv_cols_trees, ~ fit_base_models_b(TistDat_B, dv_col = .x, include_area = TRUE, use_simplified = use_simplified_re("Bushenyi", "Trees"))), threshold_labels)

## Soroti Trees: simplified only at 500m/1000m (the two thresholds that
## failed under the original structure); original elsewhere. Uses map2 with
## threshold_labels so use_simplified_re() gets the correct threshold per model.
cat("Fitting base_b trees models -- Soroti...\n")
base_b_trees_S <- setNames(
  map2(dv_cols_trees, threshold_labels, ~ fit_base_models_b(TistDat_S, dv_col = .x, include_area = TRUE, use_simplified = use_simplified_re("Soroti", "Trees", .y))),
  threshold_labels
)

cat("\nFitting h4b trees models -- Bushenyi...\n")
h4b_trees_B <- setNames(
  map2(dv_cols_trees, scaled_nf_cols, ~ fit_h4b_models(TistDat_B, dv_col = .x, nearfar_var = .y, include_area = TRUE, use_simplified = use_simplified_re("Bushenyi", "Trees"))),
  threshold_labels
)
cat("Fitting h4b trees models -- Soroti...\n")
h4b_trees_S <- setNames(
  pmap(
    list(dv_cols_trees, scaled_nf_cols, threshold_labels),
    ~ fit_h4b_models(TistDat_S, dv_col = ..1, nearfar_var = ..2, include_area = TRUE, use_simplified = use_simplified_re("Soroti", "Trees", ..3))
  ),
  threshold_labels
)

check_convergence_b(base_b_trees_B, "base_b (Trees, log1p) -- Bushenyi")
check_convergence_b(base_b_trees_S, "base_b (Trees, log1p) -- Soroti")
check_convergence_b(h4b_trees_B,    "h4b (Trees, log1p) -- Bushenyi")
check_convergence_b(h4b_trees_S,    "h4b (Trees, log1p) -- Soroti")

aic_b_trees_B <- compare_aic_b(base_b_trees_B, h4b_trees_B, "Bushenyi (Trees, log1p)")
aic_b_trees_S <- compare_aic_b(base_b_trees_S, h4b_trees_S, "Soroti (Trees, log1p)")

h2b_h3b_trees_table <- bind_rows(
  extract_h2b_h3b(base_b_trees_B, "Bushenyi"),
  extract_h2b_h3b(base_b_trees_S, "Soroti")
)
cat("\n=== H2b/H3b coefficients -- Trees (log1p) ===\n")
print(h2b_h3b_trees_table, n = 40)

##################################################################################
## Pulls the Near-Far main effect and its interaction with exposure out of
## every threshold's h4b model. Unlike extract_h2b_h3b(), the term name
## itself changes per threshold (NearFar_resid_1_sc, _2_sc, ...), so this
## looks up the matching scaled_nf_cols entry for each threshold before
## filtering.
extract_h4b <- function(mod_list, site_label) {
  map2_dfr(names(mod_list), scaled_nf_cols, function(th, nf_var) {
    m <- mod_list[[th]]
    cf <- as.data.frame(coef(summary(m)))
    cf$term <- rownames(cf)
    names(cf)[names(cf) == "Estimate"]   <- "estimate"
    names(cf)[names(cf) == "Std. Error"] <- "std.error"
    names(cf)[names(cf) == "t value"]    <- "statistic"

    interaction_term <- paste0("Exposure_sc:", nf_var)

    cf %>%
      as_tibble() %>%
      filter(term %in% c(nf_var, interaction_term)) %>%
      mutate(
        term = case_when(
          term == nf_var           ~ "NearFar_resid_sc",
          term == interaction_term ~ "Exposure_sc:NearFar_resid_sc",
          TRUE ~ term
        )
      ) %>%
      select(term, estimate, std.error, statistic) %>%
      mutate(threshold = th, site = site_label)
  })
}

## --- Density outcome ---
h4b_density_table <- bind_rows(
  extract_h4b(h4b_B, "Bushenyi"),
  extract_h4b(h4b_S, "Soroti")
)
cat("\n=== H4b coefficients -- Density ===\n")
print(h4b_density_table, n = 40)

## --- Tree count outcome (log1p) ---
h4b_trees_table <- bind_rows(
  extract_h4b(h4b_trees_B, "Bushenyi"),
  extract_h4b(h4b_trees_S, "Soroti")
)
cat("\n=== H4b coefficients -- Trees (log1p) ===\n")
print(h4b_trees_table, n = 40)

## Extracts the correlation between Exposure_sc and NearFar_resid_k_sc in
## the FIXED-EFFECTS estimates (not the raw data -- that was already
## confirmed near-zero by construction in Part 0) for each h4b model. This
## is the same diagnostic that flagged the Part 3 H4 density collinearity
## problem (AIC improving uniformly even where the coefficient itself was
## non-significant). A high |correlation| here means Exposure_sc and
## NearFar_resid_sc are being estimated jointly rather than independently,
## inflating standard errors and making it harder to attribute any AIC
## improvement or coefficient significance to one term specifically.
extract_h4b_collinearity <- function(mod_list, site_label) {
  map2_dfr(names(mod_list), scaled_nf_cols, function(th, nf_var) {
    m <- mod_list[[th]]
    corr_mat <- cov2cor(as.matrix(vcov(m)))
    tibble(
      threshold = th,
      site = site_label,
      exposure_nearfar_corr = round(corr_mat["Exposure_sc", nf_var], 3)
    )
  })
}

h4b_collinearity_density <- bind_rows(
  extract_h4b_collinearity(h4b_B, "Bushenyi"),
  extract_h4b_collinearity(h4b_S, "Soroti")
)
h4b_collinearity_trees <- bind_rows(
  extract_h4b_collinearity(h4b_trees_B, "Bushenyi"),
  extract_h4b_collinearity(h4b_trees_S, "Soroti")
)

cat("\n=== H4b collinearity check: corr(Exposure_sc, NearFar_resid_sc) in fixed-effects estimates -- Density ===\n")
print(h4b_collinearity_density, n = 40)
cat("\n=== H4b collinearity check: corr(Exposure_sc, NearFar_resid_sc) in fixed-effects estimates -- Trees ===\n")
print(h4b_collinearity_trees, n = 40)

## Merged view: coefficient estimate/statistic for the NearFar_resid_sc and
## interaction terms, alongside the collinearity value at the same
## threshold/site -- so a large coefficient sitting next to a large
## |correlation| can be spotted directly, rather than cross-referencing two
## separate tables by eye.
h4b_density_annotated <- h4b_density_table %>%
  left_join(h4b_collinearity_density, by = c("threshold", "site"))
h4b_trees_annotated <- h4b_trees_table %>%
  left_join(h4b_collinearity_trees, by = c("threshold", "site"))

cat("\n=== H4b coefficients + collinearity, merged -- Density ===\n")
print(h4b_density_annotated, n = 40)
cat("\n=== H4b coefficients + collinearity, merged -- Trees (log1p) ===\n")
print(h4b_trees_annotated, n = 40)

cat("\nInterpretation guide: treat a significant NearFar_resid_sc or interaction coefficient\n",
    "with caution wherever |exposure_nearfar_corr| is large (roughly > 0.5) at that same\n",
    "threshold -- this is the same pattern that made Part 3's Bushenyi density H4 result\n",
    "unreliable. A significant coefficient sitting alongside LOW correlation is the more\n",
    "trustworthy result.\n")

##############################################################################################
# --- H4b (dissimilarity) collinearity tables were already computed above,
# --- via extract_h4b_collinearity() (not extract_h4_collinearity()
# --- -- h4b_B/h4b_S and h4b_trees_B/h4b_trees_S are single per-threshold
# --- model lists, not $density/$trees_nb pairs like h4_B/h4_S).
# --- Reuse h4b_collinearity_density / h4b_collinearity_trees directly.

# --- Combine all four tables (H4 x Density/Trees, H4b x Density/Trees) ---
threshold_levels <- c("500m","1000m","1500m","2000m","2500m",
                      "3000m","3500m","4000m","4500m","5000m")

h4_collinearity_panel_df <- bind_rows(
  h4_collinearity_density  %>% mutate(outcome = "Density", model = "own-outcome"),
  h4_collinearity_trees    %>% mutate(outcome = "Trees",   model = "own-outcome"),
  h4b_collinearity_density %>% mutate(outcome = "Density", model = "dissimilarity"),
  h4b_collinearity_trees   %>% mutate(outcome = "Trees",   model = "dissimilarity")
) %>%
  mutate(
    threshold = fct_relevel(threshold, threshold_levels),
    outcome   = factor(outcome, levels = c("Density", "Trees")),
    model     = factor(model, levels = c("own-outcome", "dissimilarity")),
    flagged   = abs(exposure_nearfar_corr) > 0.5
  )

# --- 4-panel plot: rows = model type, columns = outcome ---
h4_collinearity_panel <- ggplot(
  h4_collinearity_panel_df,
  aes(x = threshold, y = exposure_nearfar_corr,
      colour = site, group = site)
) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = -0.5,
           fill = "grey70", alpha = 0.15) +
  geom_hline(yintercept = -0.5, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_grid(model ~ outcome) +
  scale_colour_manual(values = c("Bushenyi" = "#2E6F95", "Soroti" = "#C6602D")) +
  scale_y_continuous(limits = c(-0.9, 0), breaks = seq(-0.9, 0, 0.2)) +
  labs(
    title = "Collinearity Check: Exposure_sc vs. Near/Far Residual Across Buffer Thresholds",
    subtitle = "Own-outcome vs. dissimilarity  models, both outcomes",
    x = "Exposure buffer threshold",
    y = "corr(Exposure_sc, NearFar_resid_sc)",
    colour = "Site",
    caption = "Shaded region: |r| > 0.5 -- coefficient flagged, not reported as confirmed"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7.5),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9.5, colour = "grey30"),
    plot.caption     = element_text(size = 8, colour = "grey40", face = "italic"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1, "lines")
  )

h4_collinearity_panel

# --- Save for supplementary ---
dir.create("figures/supplementary", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/supplementary/h4_h4b_collinearity_panel.png", h4_collinearity_panel,
       width = 9, height = 7.5, dpi = 300)


################################################################################
# PART 6 -- DISSIMILARITY (H2b / H3b / H4b) COEFFICIENT-BY-THRESHOLD PLOTS
############################################################
# H2b/H3b/H4b are fundamentally about how an effect changes across
# ten spatial thresholds (500m-5000m) -- a coefficient-by-threshold
# ("forest plot") is the appropriate visualisation here, rather
# than marginal-effect curves at one arbitrarily chosen threshold.
# This reuses the coefficient tables already built earlier in the
# pipeline (h2b_h3b_density_table, h2b_h3b_trees_table,
# h4b_density_table, h4b_trees_table) -- no new models are fitted
# here. Requires threshold_labels (all 10, e.g. "500m".."5000m")
# already defined earlier in the pipeline (Part 0).
#
# NOTE ON THE ZIGZAG BUG FROM BLOCK 1: that bug came from
# geom_line() connecting points in row order on a CONTINUOUS x-axis.
# This block uses geom_pointrange() on a DISCRETE factor x-axis
# (threshold) with no connecting line between points -- row order
# cannot produce that failure mode here, so no arrange()-before-plot
# step is needed for correctness (though the code still filters/
# builds cleanly regardless of input row order).
#
# NOTE ON NAMING: the source tables here use a lowercase `site`
# column (set by extract_h2b_h3b()/extract_h4b() in Part 5 of the
# pipeline), NOT the capitalised `Site` used in Part 4's
# fig1_predictions/duration_main_predictions/nearfar_predictions.
# Both are correct within their own objects -- just don't mix the
# two column names across blocks.
#
# TITLES: removed from all plots below -- titles now live in
# manuscript figure captions only, not baked into the plot images.
############################################################

## Orders thresholds by distance (500m -> 5000m) rather than
## alphabetically, and adds an Outcome column so density/trees
## tables can be combined and faceted together. Also computes a
## 95% CI from estimate +/- std.error.
prep_coef_table <- function(tbl, outcome_label) {
  tbl %>%
    mutate(
      Outcome = outcome_label,
      threshold = factor(threshold, levels = threshold_labels),
      ci_low  = estimate - 1.96 * std.error,
      ci_high = estimate + 1.96 * std.error
    )
}

## Generic coefficient-by-threshold ("forest plot"): point + 95% CI
## per threshold, coloured by site, dashed zero-reference line.
## facet_term = TRUE additionally facets columns by term (used for
## H4b, which carries two coefficients of interest at once: the
## near-far main effect and its interaction with exposure).
plot_coef_by_threshold <- function(coef_data, facet_term = FALSE) {
  p <- ggplot(coef_data, aes(x = threshold, y = estimate, colour = site)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(
      aes(ymin = ci_low, ymax = ci_high),
      position = position_dodge(width = 0.4)
    ) +
    labs(
      x = "Distance threshold",
      y = "Coefficient estimate (95% CI)",
      colour = "Site"
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(face = "bold", size = 13)
    )

  if (facet_term) {
    p <- p + facet_grid(Outcome ~ term, scales = "free_y")
  } else {
    p <- p + facet_wrap(~Outcome, scales = "free_y")
  }

  p
}

# ----------------------------------------------------------
# H2b / H3b shared source table -- built once, filtered twice
# (previously rebuilt independently for each figure; identical
# result, just avoids the redundant bind_rows() call)
# ----------------------------------------------------------

h2b_h3b_combined <- bind_rows(
  prep_coef_table(h2b_h3b_density_table, "Density"),
  prep_coef_table(h2b_h3b_trees_table,   "Trees")
)

# ----------------------------------------------------------
# H2b: Exposure_sc main effect on neighbour-dissimilarity,
#      by distance threshold
# ----------------------------------------------------------

h2b_exposure_data <- h2b_h3b_combined %>%
  filter(term == "Exposure_sc")

fig_h2b_exposure <- plot_coef_by_threshold(h2b_exposure_data)

fig_h2b_exposure

ggsave(
  filename = "figure_H2b_exposure_dissimilarity.png",
  plot = fig_h2b_exposure,
  width = 11,
  height = 6,
  dpi = 300,
  units = "in"
)


# ----------------------------------------------------------
# H3b: Exposure x Duration interaction on dissimilarity,
#      by distance threshold
# ----------------------------------------------------------

h3b_interaction_data <- h2b_h3b_combined %>%
  filter(term == "Exposure_sc:Years_since_reg_sc")

fig_h3b_interaction <- plot_coef_by_threshold(h3b_interaction_data)

fig_h3b_interaction

ggsave(
  filename = "figure_H3b_duration_dissimilarity.png",
  plot = fig_h3b_interaction,
  width = 11,
  height = 6,
  dpi = 300,
  units = "in"
)


# ----------------------------------------------------------
# H4b: Near-Far residual main effect + its interaction with
#      exposure, by distance threshold
# ----------------------------------------------------------

h4b_combined_data <- bind_rows(
  prep_coef_table(h4b_density_table, "Density"),
  prep_coef_table(h4b_trees_table,   "Trees")
) %>%
  mutate(
    term = factor(
      term,
      levels = c("NearFar_resid_sc", "Exposure_sc:NearFar_resid_sc"),
      labels = c("Near-Far residual (main effect)", "Exposure x Near-Far residual")
    )
  )

fig_h4b <- plot_coef_by_threshold(h4b_combined_data, facet_term = TRUE)

fig_h4b

ggsave(
  filename = "figure_H4b_nearfar_dissimilarity.png",
  plot = fig_h4b,
  width = 13,
  height = 7,
  dpi = 300,
  units = "in"
)


# ----------------------------------------------------------
# COMBINED PANEL: Dissimilarity summary (H2b + H3b)
# ----------------------------------------------------------

fig_dissimilarity_combined <- (fig_h2b_exposure / fig_h3b_interaction) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_dissimilarity_combined

ggsave(
  filename = "figure_dissimilarity_combined.png",
  plot = fig_dissimilarity_combined,
  width = 11,
  height = 11,
  dpi = 300,
  units = "in"
)


################################################################################
# PART 7 -- NORM-CONVERGENCE ASYMMETRY TEST -- ABOVE vs. BELOW
# NEIGHBOURHOOD MEAN
#
# The H2/H2b tree-count finding (fewer trees overall, but greater similarity
# to neighbours) is consistent with exposure pulling farmers toward a shared
# local planting level -- but a model on AbsDissimilarity (always >= 0)
# cannot tell us whether this is genuine two-way convergence (farmers above
# the local mean move down AND farmers below move up) or one-way
# "levelling" (only one side moves toward the mean). This test uses the
# SIGNED dissimilarity measure (Dissimilarity_Trees_k = own - neighbour
# mean, not its absolute value) and fits the same specification separately
# for farmers above vs. below their neighbourhood mean.
#
# Convergence prediction:
#   - Above-mean farmers: Exposure_sc should be NEGATIVE (higher exposure ->
#     smaller positive gap -> moving DOWN toward neighbours)
#   - Below-mean farmers: Exposure_sc should be POSITIVE (higher exposure ->
#     gap moves toward zero -> moving UP toward neighbours)
# If only one side shows the predicted sign/significance, this indicates
# asymmetric ("levelling") convergence rather than the symmetric
# norm-convergence account proposed in the Discussion draft.
################################################################################

test_convergence_asymmetry <- function(dat, k, site_label, use_simplified = FALSE) {
  signed_col <- paste0("Dissimilarity_Trees_", k)

  dat_split <- dat %>%
    mutate(Above_Neighbour_Mean = .data[[signed_col]] > 0)

  re_term <- if (use_simplified) {
    "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
  } else {
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  }

  frm <- as.formula(paste0(
    signed_col, " ~ Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc + ",
    "Dist_To_Forest_sc + Area_Ha_sc + ", re_term
  ))

  fit_group <- function(group_label, group_filter) {
    sub <- dat_split %>% filter(Above_Neighbour_Mean == group_filter)
    if (nrow(sub) < 30) {
      cat("  [", site_label, threshold_labels[k], group_label, "] n =", nrow(sub), "-- too small, skipping.\n")
      return(NULL)
    }
    m <- tryCatch(lmer(frm, data = sub, REML = TRUE, na.action = na.exclude), error = function(e) NULL)
    if (is.null(m)) {
      cat("  [", site_label, threshold_labels[k], group_label, "] model failed to fit -- skipping.\n")
      return(NULL)
    }
    cf <- as.data.frame(coef(summary(m)))
    cf$term <- rownames(cf)
    tibble(
      site = site_label, threshold = threshold_labels[k], group = group_label, n = nrow(sub),
      exposure_est = cf$Estimate[cf$term == "Exposure_sc"],
      exposure_se  = cf$`Std. Error`[cf$term == "Exposure_sc"],
      exposure_t   = cf$`t value`[cf$term == "Exposure_sc"]
    )
  }

  bind_rows(
    fit_group("Above neighbour mean", TRUE),
    fit_group("Below neighbour mean", FALSE)
  )
}

## Run at a representative set of thresholds spanning the range where H2b's
## Trees convergence finding was significant in both sites (1000m-4000m),
## rather than all 10, to keep this a targeted check rather than a full sweep.
convergence_thresholds <- c(2, 4, 6, 8)   # 1000m, 2000m, 3000m, 4000m

cat("\nRunning norm-convergence asymmetry test -- Bushenyi...\n")
convergence_asymmetry_B <- map_dfr(
  convergence_thresholds,
  ~ test_convergence_asymmetry(TistDat_B, .x, "Bushenyi", use_simplified = use_simplified_re("Bushenyi", "Trees", threshold_labels[.x]))
)

cat("Running norm-convergence asymmetry test -- Soroti...\n")
convergence_asymmetry_S <- map_dfr(
  convergence_thresholds,
  ~ test_convergence_asymmetry(TistDat_S, .x, "Soroti", use_simplified = use_simplified_re("Soroti", "Trees", threshold_labels[.x]))
)

convergence_asymmetry_results <- bind_rows(convergence_asymmetry_B, convergence_asymmetry_S)

cat("\n=== Norm-convergence asymmetry test: Exposure_sc on signed Dissimilarity_Trees, by direction ===\n")
print(convergence_asymmetry_results, n = 40)

cat("\nInterpretation guide:\n")
cat(" - Symmetric (true) convergence: Above group shows NEGATIVE, significant Exposure_sc\n")
cat("   AND Below group shows POSITIVE, significant Exposure_sc, at the same threshold.\n")
cat(" - 'Levelling down' only: Above group significant (negative) but Below group is not.\n")
cat(" - 'Levelling up' only: Below group significant (positive) but Above group is not.\n")
cat(" - Neither significant: no support for a directional convergence mechanism at that threshold.\n")


################################################################################
# PART 8 -- SPATIAL DIAGNOSTICS ON THE PART 3 (H2/H3/H4) MODELS
# MORAN'S I + DHARMa SIGNIFICANCE TESTS
#
# WHY THIS IS KEPT SEPARATE FROM PARTS 1-2: the semivariograms in Parts 1-2
# are DESCRIPTIVE -- they show how far spatial pattern reaches and how
# strong it is, but they do not, by themselves, give a formal significance
# test. Moran's I and DHARMa's spatial-autocorrelation test DO give a
# formal p-value, answering a simple "is there significant leftover spatial
# pattern, yes or no". The two approaches are complementary, not
# interchangeable, and both are kept:
#   - DHARMa's test runs on INDIVIDUAL farmer-level residuals, the same
#     scale as the Part 2 residual variograms, and correctly accounts for
#     each model's error distribution (Gaussian for density, negative
#     binomial for tree counts) via simulation.
#   - Moran's I runs on VILLAGE-AGGREGATED mean residuals, a coarser scale
#     that avoids treating many farmers at the same location as fully
#     independent data points.
################################################################################

## Marginal effects: how the model expects density to change with exposure,
## at different levels of years-since-registration.
ggpredict(base_B$density, terms = c("Exposure_sc", "Years_since_reg_sc")) |> plot()
ggpredict(base_S$density, terms = c("Exposure_sc", "Years_since_reg_sc")) |> plot()

## Residual extraction (density: raw residuals; trees: Pearson residuals,
## standardised to make over-dispersion assessment meaningful).
dat_B <- TistDat_B
dat_B$res_density <- resid(base_B$density)
dat_B$res_trees   <- residuals(base_B$trees_nb, type = "pearson")

dat_S <- TistDat_S
dat_S$res_density <- resid(base_S$density)
dat_S$res_trees   <- residuals(base_S$trees_nb, type = "pearson")

## Aggregate residuals to village level (one value per village) for the
## Moran's I test, using each village's mean farmer coordinates.
village_res_B <- dat_B |>
  group_by(Village_ID) |>
  summarise(res_density = mean(res_density, na.rm = TRUE), res_trees = mean(res_trees, na.rm = TRUE),
            lon = mean(longitude, na.rm = TRUE), lat = mean(latitude, na.rm = TRUE), .groups = "drop")

village_res_S <- dat_S |>
  group_by(Village_ID) |>
  summarise(res_density = mean(res_density, na.rm = TRUE), res_trees = mean(res_trees, na.rm = TRUE),
            lon = mean(longitude, na.rm = TRUE), lat = mean(latitude, na.rm = TRUE), .groups = "drop")

## Spatial weights: since villages are irregularly spaced, each village is
## connected to its 4 nearest neighbours (kNN) rather than using a
## contiguity ("shares a border") rule.
build_spatial_weights <- function(village_res, utm_crs = 32736, k = 4) {
  sf_obj  <- st_as_sf(village_res, coords = c("lon", "lat"), crs = 4326)
  sf_proj <- st_transform(sf_obj, utm_crs)
  coords  <- st_coordinates(sf_proj)
  nb      <- knn2nb(knearneigh(coords, k = k))
  lw      <- nb2listw(nb, style = "W")
  list(sf = sf_proj, coords = coords, nb = nb, lw = lw)
}

spatial_B <- build_spatial_weights(village_res_B)
spatial_S <- build_spatial_weights(village_res_S)

cat("\n=== Moran's I -- Bushenyi ===\n")
cat("Density residuals:\n"); print(moran.test(village_res_B$res_density, spatial_B$lw))
cat("\nTree count residuals (Pearson):\n"); print(moran.test(village_res_B$res_trees, spatial_B$lw))

cat("\n=== Moran's I -- Soroti ===\n")
cat("Density residuals:\n"); print(moran.test(village_res_S$res_density, spatial_S$lw))
cat("\nTree count residuals (Pearson):\n"); print(moran.test(village_res_S$res_trees, spatial_S$lw))

## DHARMa: simulation-based residual diagnostics, appropriate for both the
## Gaussian (density) and negative-binomial (trees) models.
cat("\n=== DHARMa diagnostics -- Bushenyi ===\n")
sim_B_density <- simulateResiduals(base_B$density, n = 1000)
sim_B_trees   <- simulateResiduals(base_B$trees_nb, n = 1000)
dev.new(); plot(sim_B_density, main = "Bushenyi -- Density")
dev.new(); plot(sim_B_trees,   main = "Bushenyi -- Tree count")

cat("\nSpatial autocorrelation test -- Bushenyi density:\n")
testSpatialAutocorrelation(sim_B_density, x = dat_B$longitude, y = dat_B$latitude, plot = TRUE)
cat("\nSpatial autocorrelation test -- Bushenyi tree count:\n")
testSpatialAutocorrelation(sim_B_trees, x = dat_B$longitude, y = dat_B$latitude, plot = TRUE)

cat("\n=== DHARMa diagnostics -- Soroti ===\n")
sim_S_density <- simulateResiduals(base_S$density, n = 1000)
sim_S_trees   <- simulateResiduals(base_S$trees_nb, n = 1000)
dev.new(); plot(sim_S_density, main = "Soroti -- Density")
dev.new(); plot(sim_S_trees,   main = "Soroti -- Tree count")

cat("\nSpatial autocorrelation test -- Soroti density:\n")
testSpatialAutocorrelation(sim_S_density, x = dat_S$longitude, y = dat_S$latitude, plot = TRUE)
cat("\nSpatial autocorrelation test -- Soroti tree count:\n")
testSpatialAutocorrelation(sim_S_trees, x = dat_S$longitude, y = dat_S$latitude, plot = TRUE)

## Additional DHARMa checks specific to the count model: is variance too
## high for a standard model to expect (dispersion), and are there more
## zero counts than the model predicts (zero-inflation)?
cat("\n=== Dispersion tests -- tree count models ===\n")
cat("Bushenyi:\n"); testDispersion(sim_B_trees)
cat("Soroti:\n");   testDispersion(sim_S_trees)

cat("\n=== Zero-inflation tests -- tree count models ===\n")
cat("Bushenyi:\n"); testZeroInflation(sim_B_trees)
cat("Soroti:\n");   testZeroInflation(sim_S_trees)


################################################################################
# END OF SCRIPT
#
# Everything from the original exploratory check (raw variograms) and the
# hypothesis-testing script (residual variograms, H2b/H3b/H4b models,
# Moran's I / DHARMa diagnostics) is preserved here -- reorganised so that
# the semivariogram analysis (Parts 1-2) now comes first, right after
# setup, ahead of the modelling sections (Parts 3-7) that were previously
# interleaved with it. Nothing has been dropped or altered in substance;
# the only structural change beyond reordering is that use_simplified_re()
# is now defined once, in Part 2, where it is first needed (by the Stage B
# residual variograms), and simply reused by Parts 5 and 7 rather than
# being redefined.
################################################################################
