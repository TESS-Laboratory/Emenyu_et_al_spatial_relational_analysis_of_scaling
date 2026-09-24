################################################################################
# TIST MODELLING PIPELINE -- Chapter 3
#
#   PART 0 -- Setup: paths, settings, helpers, data, scaling, site split
#   PART 1 -- Semivariograms, Stage A: raw outcomes
#   PART 2 -- Semivariograms, Stage B: residuals after exposure alone
#   PART 3 -- Main models (H2/H3/H4) and model-fit outputs
#   PART 4 -- Diagnostics and figures for the main models
#   PART 5 -- Neighbour-dissimilarity models (H2b/H3b/H4b)
#   PART 6 -- Coefficient-by-threshold figures
#   PART 7 -- Norm-convergence asymmetry test
#   PART 8 -- Residual diagnostics: Moran's I, DHARMa, zero-truncation check
#   PART 9 -- RESULTS SUMMARY: every result printed together and saved
#
# Every result is stored with store_result() as it is produced. Part 9
# prints them all, in order, and writes the same output to
# <results_dir>/results_summary.txt. Tables are also saved as CSV.
#
# CHANGES FROM THE PREVIOUS VERSION
#   - One set of output paths, defined once (previously defined twice).
#   - The dataset is read once; the variogram sample is derived from it
#     (previously the same CSV was read a second time).
#   - Model formulas are built by one helper (build_formula), so the base,
#     H4, dissimilarity, convergence and zero-truncated models cannot drift
#     apart.
#   - One coefficient extractor (tidy_fixed) and one collinearity extractor
#     replace the four near-identical versions.
#   - DHARMa tests run without on-screen plots (the cause of the earlier
#     stalls); plots can be saved to PNG by setting save_dharma_plots.
#   - Added: zero-truncated negative binomial robustness check (Part 8).
#   - Removed: the ggpredict() marginal-effect plots at the start of Part 8
#     (duplicated Figure 3.10) and the uncorrected tenure-tercile
#     predictions (tenure-tercile predictions are now reported once, for a typical group).
#   - The magick image objects are no longer named `c`, which masked base c().
################################################################################

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


################################################################################
# PART 0 -- SETUP
################################################################################

## ---- Paths (defined once) ---------------------------------------------------
data_path     <- "Data/TISTDat/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv"
output_dir    <- "Output/Manuscript 3 graphs"
variogram_dir <- file.path(output_dir, "variograms")
fig_dir       <- file.path(output_dir, "figures")
supp_fig_dir  <- file.path(output_dir, "figures", "supplementary")
results_dir   <- file.path(output_dir, "model_results")
for (d in c(variogram_dir, fig_dir, supp_fig_dir, results_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- Settings ----------------------------------------------------------------
## Variogram sample. TRUE reproduces the published Table 3.2 / Figure 3.7,
## which dropped every row with a missing value in ANY column (36,565
## records). FALSE drops rows only where coordinates or the outcome are
## missing, matching the modelled sample more closely.
variogram_complete_cases <- TRUE

## DHARMa diagnostic plots are slow with ~3,000 points and 1,000
## simulations and can stall RStudio's plot pane. FALSE = tests only.
## TRUE = plots are written to PNG files (never to the screen).
save_dharma_plots <- FALSE
n_dharma_sim      <- 1000

## DHARMa is simulation-based; fixing the seed makes its p-values
## reproducible. (Values may differ very slightly from earlier unseeded runs.)
set.seed(2026)

options(pillar.sigfig = 5)

## ---- Results store -----------------------------------------------------------
## Every result is added here in the order produced, then printed in Part 9.
results <- list()
store_result <- function(name, obj, csv = NULL) {
  if (is.data.frame(obj)) obj <- as_tibble(obj)
  results[[name]] <<- obj
  if (!is.null(csv) && is.data.frame(obj)) {
    write.csv(obj, file.path(results_dir, csv), row.names = FALSE)
  }
  invisible(obj)
}

## ---- Shared helpers ----------------------------------------------------------

## Tidy fixed effects for lmer (density, dissimilarity) and glmmTMB (tree
## count) models.
##   p.value: glmmTMB = Wald z-test; lmer = normal approximation to t
##            (the same basis as the |t| > 2 rule in Table S8).
##   pct_change: 100 * (exp(estimate) - 1); meaningful for log-link
##               tree-count models only (NA for lmer).
tidy_fixed <- function(m, model_label) {
  if (inherits(m, "glmmTMB")) {
    cf <- summary(m)$coefficients$cond
    tibble(
      model      = model_label,
      term       = rownames(cf),
      estimate   = cf[, "Estimate"],
      std.error  = cf[, "Std. Error"],
      statistic  = cf[, "z value"],
      p.value    = cf[, "Pr(>|z|)"],
      pct_change = 100 * (exp(cf[, "Estimate"]) - 1)
    )
  } else {
    cf <- coef(summary(m))
    tibble(
      model      = model_label,
      term       = rownames(cf),
      estimate   = cf[, "Estimate"],
      std.error  = cf[, "Std. Error"],
      statistic  = cf[, "t value"],
      p.value    = 2 * pnorm(-abs(cf[, "t value"])),
      pct_change = NA_real_
    )
  }
}

## One row per hypothesis test (Moran's I, DHARMa, chi-square).
htest_row <- function(h, label) {
  vals <- c(h$statistic, h$estimate)
  tibble(
    test       = label,
    method     = h$method,
    statistics = paste(names(vals), signif(vals, 4), sep = " = ", collapse = "; "),
    p.value    = h$p.value
  )
}

## Variance-covariance matrix of the fixed effects, for lmer or glmmTMB.
get_vcov_matrix <- function(m) {
  v <- vcov(m)
  if (is.list(v) && !is.matrix(v)) v <- v$cond
  as.matrix(v)
}

## ---- Model formulas (one builder for every model in the pipeline) -----------
re_full   <- "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
re_simple <- "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
core_fixed <- "Exposure_sc + Years_since_reg_sc + Exposure_sc:Years_since_reg_sc + Dist_To_Forest_sc"

build_formula <- function(dv, area = FALSE, nearfar = NULL, simplified = FALSE) {
  as.formula(paste0(
    dv, " ~ ", core_fixed,
    if (area) " + Area_Ha_sc" else "",
    if (!is.null(nearfar)) paste0(" + ", nearfar, " + Exposure_sc:", nearfar) else "",
    " + ", if (simplified) re_simple else re_full
  ))
}

## Which random-effects structure to use, by site / outcome / threshold.
## The full structure produced singular fits for the density models and for
## Soroti tree-count dissimilarity at 500 m and 1,000 m; dropping
## Admin_Districts/Subcounty resolves these at negligible AIC cost. For
## Bushenyi tree count the full structure never failed, and simplifying it
## worsens AIC (~+147), so it is never simplified there. Applies to the
## Part 2 exposure-only models and the Part 5 and Part 7 models; the Part 3
## main models always use the full structure (Table S8).
use_simplified_re <- function(site, outcome, threshold_label = NA) {
  if (outcome == "Density") return(TRUE)
  if (outcome == "Trees" && site == "Bushenyi") return(FALSE)
  if (outcome == "Trees" && site == "Soroti") return(threshold_label %in% c("500m", "1000m"))
  FALSE
}

## ---- Load, scale, split ------------------------------------------------------
TistDat_H <- read.csv(data_path) %>%
  distinct() %>%
  mutate(
    Dist_To_Forest_sc  = as.numeric(scale(Dist_To_Forest_m)),
    Area_Ha_sc         = as.numeric(scale(Area_Ha)),
    Trees_sc           = as.numeric(scale(Trees)),
    Exposure_sc        = as.numeric(scale(Exposure)),
    Years_since_reg_sc = as.numeric(scale(Years_since_reg)),
    ## Spatial-concentration (Near-Far) indices at 10 thresholds, 500-5,000 m.
    across(NearFar_resid_1:NearFar_resid_10, ~ as.numeric(scale(.x)), .names = "{.col}_sc")
  )
## NOTE: predictors are standardised on the combined two-site dataset.

scaled_nf_cols   <- paste0("NearFar_resid_", 1:10, "_sc")
threshold_labels <- paste0(seq(500, 5000, by = 500), "m")

store_result("0.1 Near-Far indices: correlation with Exposure_sc (should be ~0)",
             tibble(threshold = threshold_labels,
                    corr_with_exposure = sapply(scaled_nf_cols, function(v)
                      round(cor(TistDat_H[[v]], TistDat_H$Exposure_sc, use = "complete.obs"), 8))))

## Every model is fitted separately per site.
TistDat_B <- TistDat_H |> filter(Proj_Area == "Bushenyi")
TistDat_S <- TistDat_H |> filter(Proj_Area == "Soroti")

## Tenure terciles, cut per site. Cut-points fall on whole registration
## years, so the three groups are not equal in size.
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

## ---- Sample sizes and summary statistics -------------------------------------
store_result("0.2 Sample sizes (cleaned dataset and model samples)", tibble(
  sample          = c("Cleaned dataset (all Uganda)", "Bushenyi model sample", "Soroti model sample"),
  records         = c(nrow(TistDat_H), nrow(TistDat_B), nrow(TistDat_S)),
  unique_farmers  = c(n_distinct(TistDat_H$Farmer_ID), n_distinct(TistDat_B$Farmer_ID), n_distinct(TistDat_S$Farmer_ID)),
  unique_groups   = c(n_distinct(TistDat_H$Group_ID),   n_distinct(TistDat_B$Group_ID),   n_distinct(TistDat_S$Group_ID)),
  unique_villages = c(n_distinct(TistDat_H$Village_ID), n_distinct(TistDat_B$Village_ID), n_distinct(TistDat_S$Village_ID))
), csv = "sample_sizes.csv")

store_result("0.3 Tenure tercile counts", bind_rows(
  count(TistDat_B, Duration_Tercile) %>% mutate(site = "Bushenyi"),
  count(TistDat_S, Duration_Tercile) %>% mutate(site = "Soroti")
) %>% select(site, everything()))

summarise_model_sample <- function(dat, site_label) {
  vars <- c("Dist_To_Forest_m", "Exposure", "Density_winsor99", "Trees", "Area_Ha", "Years_since_reg")
  map_dfr(vars, function(v) {
    x <- dat[[v]]
    tibble(
      site = site_label, variable = v, n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
      p90 = unname(quantile(x, 0.9, na.rm = TRUE)),
      min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)
    )
  })
}
store_result("0.4 Model-sample summary statistics (Tables S6-S7)",
             bind_rows(summarise_model_sample(TistDat_B, "Bushenyi"),
                       summarise_model_sample(TistDat_S, "Soroti")),
             csv = "model_sample_summary.csv")


################################################################################
# PART 1 -- SEMIVARIOGRAMS, STAGE A: RAW OUTCOMES
#
# Do farms close together resemble each other more than farms far apart?
# Semivariograms fitted to the raw (log1p) outcomes, per site x outcome.
################################################################################

## Fit a semivariogram to any numeric column (raw outcome or residuals).
## Retries once on a degenerate fit, then falls back to a pure-nugget model.
fit_variogram_generic <- function(data_with_coords, value_col, cutoff = 8000, width = 500) {
  data_clean <- data_with_coords %>% filter(!is.na(.data[[value_col]]))
  data_sf <- st_as_sf(data_clean, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(32636)   # UTM 36N, distances in metres
  data_sp <- as(data_sf, "Spatial")
  
  vgm_empirical <- variogram(as.formula(paste0(value_col, " ~ 1")), data_sp, cutoff = cutoff, width = width)
  
  fit_attempt <- function(start_range) {
    tryCatch(
      fit.variogram(vgm_empirical,
                    vgm(psill = var(data_sp[[value_col]], na.rm = TRUE), model = "Sph",
                        range = start_range, nugget = 0)),
      error = function(e) NULL
    )
  }
  range_ok <- function(v) !is.null(v) && any(v$model == "Sph") &&
    is.finite(v$range[v$model == "Sph"]) && v$range[v$model == "Sph"] > 0
  
  vgm_fitted <- fit_attempt(2000)
  if (!range_ok(vgm_fitted)) {
    message("  Spherical fit failed for '", value_col, "' -- retrying with a different starting range.")
    vgm_fitted <- fit_attempt(max(vgm_empirical$dist) / 2)
  }
  if (!range_ok(vgm_fitted)) {
    message("  Spherical fit still failed for '", value_col,
            "' -- pure-nugget fallback. Treat 0% spatial dependence as 'not estimable'.")
    vgm_fitted <- vgm(psill = var(data_sp[[value_col]], na.rm = TRUE), model = "Nug")
  }
  
  nugget_val       <- if (any(vgm_fitted$model == "Nug")) vgm_fitted$psill[vgm_fitted$model == "Nug"] else 0
  partial_sill_val <- if (any(vgm_fitted$model == "Sph")) vgm_fitted$psill[vgm_fitted$model == "Sph"] else 0
  range_val        <- if (any(vgm_fitted$model == "Sph")) vgm_fitted$range[vgm_fitted$model == "Sph"] else 0
  total_sill_val   <- nugget_val + partial_sill_val
  
  list(
    empirical = vgm_empirical, fitted = vgm_fitted,
    fit_curve = variogramLine(vgm_fitted, maxdist = max(vgm_empirical$dist)),
    nugget = nugget_val, partial_sill = partial_sill_val,
    total_sill = total_sill_val, range_m = range_val,
    spatial_dependence = if (total_sill_val > 0) partial_sill_val / total_sill_val else 0
  )
}

## Draw one semivariogram panel on the open device.
draw_variogram_panel <- function(vgm_result, plot_title, label = NULL, show_stats = FALSE) {
  dist_full  <- c(0, vgm_result$empirical$dist)
  gamma_full <- c(vgm_result$nugget, vgm_result$empirical$gamma)
  y_max      <- max(gamma_full, vgm_result$total_sill * 1.1)
  
  plot(dist_full, gamma_full, type = "b", pch = 19, col = "black",
       xlim = c(0, max(dist_full)), ylim = c(0, y_max),
       xlab = "Distance (m)", ylab = "Semivariance", main = plot_title)
  lines(vgm_result$fit_curve$dist, vgm_result$fit_curve$gamma, col = "gray40", lwd = 2)
  abline(v = vgm_result$range_m, lty = 2, lwd = 2, col = "red")
  text(x = vgm_result$range_m, y = y_max * 0.85,
       labels = paste0("Range = ", round(vgm_result$range_m, 0), " m"),
       col = "red", cex = 1.1, pos = 4)
  
  if (show_stats) {
    spatial_pct <- round(vgm_result$spatial_dependence * 100, 0)
    text(x = par("usr")[2] - 0.03 * diff(par("usr")[1:2]), y = y_max * 0.15,
         labels = paste0("Spatial dependence = ", spatial_pct, "%\n(", 100 - spatial_pct, "% unexplained)"),
         col = "black", cex = 1.0, adj = c(1, 0))
  }
  if (!is.null(label)) {
    text(x = par("usr")[1], y = par("usr")[4], labels = label, font = 2, cex = 1.6,
         adj = c(-0.3, 1.3), xpd = NA)
  }
}

plot_variogram_standard <- function(vgm_result, plot_title, save_path, show_stats = FALSE) {
  png(save_path, width = 2400, height = 1800, res = 300)
  par(mar = c(5, 5, 3, 2), cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5)
  draw_variogram_panel(vgm_result, plot_title, show_stats = show_stats)
  dev.off()
}

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

run_raw_variogram <- function(full_data, site_name, outcome_col, outcome_label,
                              log_transform = TRUE, show_stats = FALSE) {
  site_data <- full_data %>% filter(Proj_Area == site_name)
  value_col <- outcome_col
  if (log_transform) {
    value_col <- paste0("log_", outcome_col)
    site_data[[value_col]] <- log1p(site_data[[outcome_col]])
  }
  vgm_result <- fit_variogram_generic(site_data, value_col)
  plot_variogram_standard(vgm_result, paste0(site_name, " - ", outcome_label),
                          file.path(variogram_dir, paste0(site_name, "_variogram_", outcome_col, ".png")),
                          show_stats = show_stats)
  vgm_result
}

## Variogram sample, derived from the dataset already loaded (see setting in Part 0).
TistDist_Dat <- if (variogram_complete_cases) {
  TistDat_H %>% drop_na()
} else {
  TistDat_H %>% drop_na(longitude, latitude, Density_winsor99, Trees)
}
store_result("1.1 Variogram sample size", tibble(
  complete_cases_rule = variogram_complete_cases,
  records = nrow(TistDist_Dat),
  bushenyi = sum(TistDist_Dat$Proj_Area == "Bushenyi"),
  soroti   = sum(TistDist_Dat$Proj_Area == "Soroti")
))

raw_density_Bushenyi <- run_raw_variogram(TistDist_Dat, "Bushenyi", "Density_winsor99", "Tree density")
raw_density_Soroti   <- run_raw_variogram(TistDist_Dat, "Soroti",   "Density_winsor99", "Tree density")
raw_trees_Bushenyi   <- run_raw_variogram(TistDist_Dat, "Bushenyi", "Trees",            "Tree count")
raw_trees_Soroti     <- run_raw_variogram(TistDist_Dat, "Soroti",   "Trees",            "Tree count")

raw_vgms <- list(raw_density_Bushenyi, raw_density_Soroti, raw_trees_Bushenyi, raw_trees_Soroti)

make_variogram_panel_figure(
  results_list = raw_vgms,
  titles       = c("Bushenyi - Tree density", "Soroti - Tree density",
                   "Bushenyi - Tree count",  "Soroti - Tree count"),
  save_path    = file.path(variogram_dir, "raw_variograms_panel.png"),
  show_stats   = TRUE
)

summarise_vgms <- function(vgm_list) {
  tibble(
    site       = c("Bushenyi", "Soroti", "Bushenyi", "Soroti"),
    outcome    = c("Density (log)", "Density (log)", "Trees (log)", "Trees (log)"),
    nugget     = map_dbl(vgm_list, "nugget"),
    total_sill = map_dbl(vgm_list, "total_sill"),
    range_m    = map_dbl(vgm_list, "range_m"),
    spatial_dependence_pct = round(100 * map_dbl(vgm_list, "spatial_dependence"), 0)
  )
}
raw_variogram_summary <- summarise_vgms(raw_vgms)
store_result("1.2 Raw semivariogram parameters (Table 3.2)", raw_variogram_summary,
             csv = "raw_variogram_summary.csv")

## Labelled a/b/c/d composite of the four single-plot PNGs (Figure 3.7).
panel_imgs <- map2(
  c("Bushenyi_variogram_Density_winsor99.png", "Soroti_variogram_Density_winsor99.png",
    "Bushenyi_variogram_Trees.png", "Soroti_variogram_Trees.png"),
  c("a", "b", "c", "d"),
  ~ image_annotate(image_read(file.path(variogram_dir, .x)), .y,
                   size = 80, location = "+20+20", weight = 700)
)
combined_variogram_panels <- image_append(c(
  image_append(c(panel_imgs[[1]], panel_imgs[[2]])),
  image_append(c(panel_imgs[[3]], panel_imgs[[4]]))
), stack = TRUE)
image_write(combined_variogram_panels, file.path(variogram_dir, "combined_variogram_panels.png"))


################################################################################
# PART 2 -- SEMIVARIOGRAMS, STAGE B: RESIDUALS AFTER EXPOSURE ALONE
#
# How much of the raw spatial clustering does exposure account for on its
# own? The exposure-only model is fitted on the same log1p scale as Stage A,
# so the sills are comparable: (sill_raw - sill_residual) / sill_raw.
################################################################################

run_exposure_residual_variogram <- function(dat, outcome_col, use_simplified = FALSE) {
  dat_complete <- dat %>%
    filter(!is.na(longitude), !is.na(latitude), !is.na(.data[[outcome_col]]))
  log_col <- paste0("log_", outcome_col)
  dat_complete[[log_col]] <- log1p(dat_complete[[outcome_col]])
  
  frm <- as.formula(paste0(log_col, " ~ Exposure_sc + ", if (use_simplified) re_simple else re_full))
  m_exposure <- lmer(frm, data = dat_complete, REML = TRUE, na.action = na.exclude)
  
  dat_res <- dat_complete %>% mutate(resid_val = resid(m_exposure))
  fit_variogram_generic(dat_res, "resid_val")
}

## Threshold only matters for Soroti tree count (see use_simplified_re()).
exposure_stage_b_threshold <- "500m"

exp_resid_vgms <- list(
  run_exposure_residual_variogram(TistDat_B, "Density_winsor99",
                                  use_simplified_re("Bushenyi", "Density", exposure_stage_b_threshold)),
  run_exposure_residual_variogram(TistDat_S, "Density_winsor99",
                                  use_simplified_re("Soroti", "Density", exposure_stage_b_threshold)),
  run_exposure_residual_variogram(TistDat_B, "Trees",
                                  use_simplified_re("Bushenyi", "Trees", exposure_stage_b_threshold)),
  run_exposure_residual_variogram(TistDat_S, "Trees",
                                  use_simplified_re("Soroti", "Trees", exposure_stage_b_threshold))
)

make_variogram_panel_figure(
  results_list = exp_resid_vgms,
  titles       = c("Bushenyi - Tree density (residual)", "Soroti - Tree density (residual)",
                   "Bushenyi - Tree count (residual)",  "Soroti - Tree count (residual)"),
  save_path    = file.path(variogram_dir, "exposure_residual_variograms_panel.png"),
  show_stats   = TRUE
)

exposure_residual_summary <- summarise_vgms(exp_resid_vgms)
store_result("2.1 Exposure-only residual semivariogram parameters", exposure_residual_summary,
             csv = "exposure_residual_variogram_summary.csv")

variance_explained_by_exposure <- raw_variogram_summary %>%
  select(site, outcome, sill_raw = total_sill) %>%
  mutate(sill_residual = exposure_residual_summary$total_sill,
         pct_variance_explained_by_exposure = round(100 * (sill_raw - sill_residual) / sill_raw, 1))
store_result("2.2 Variance explained by exposure alone (Section 3.9.2)", variance_explained_by_exposure,
             csv = "variance_explained_by_exposure.csv")

## Word table. Rows where the residual range did not plateau within the
## 8,000 m cutoff are shown as "Not resolved" -- set by inspecting the
## Stage B panel figure; update if the data change.
range_resolved <- c(TRUE, FALSE, FALSE, FALSE)

residual_comparison_table <- raw_variogram_summary %>%
  transmute(
    Site = site,
    Outcome = str_remove(outcome, " \\(log\\)"),
    `Range, raw` = paste0(round(range_m, 0), " m"),
    `Range, residual` = ifelse(range_resolved,
                               paste0(round(exposure_residual_summary$range_m, 0), " m"),
                               "Not resolved\u2020"),
    `Spatial dependence, raw (%)` = spatial_dependence_pct,
    `Spatial dependence, residual (%)` = exposure_residual_summary$spatial_dependence_pct,
    `Variance explained by exposure (%)` = variance_explained_by_exposure$pct_variance_explained_by_exposure
  )

ft <- flextable(residual_comparison_table) %>%
  set_caption("Comparison of raw and exposure-only residual semivariogram parameters, by site and outcome.") %>%
  add_footer_lines(
    "\u2020 Fitted spherical range exceeded the maximum tested distance (8,000 m); the empirical semivariogram had not plateaued within the observed range. Spatial dependence values for these rows should be interpreted with caution; variance explained by exposure (final column) does not share this limitation, as it is computed from total sills rather than the plateau location."
  ) %>%
  theme_vanilla() %>% autofit() %>% fontsize(size = 9, part = "all") %>% bold(part = "header")
save_as_docx(ft, path = file.path(results_dir, "Table_residual_variogram_comparison.docx"))


################################################################################
# PART 3 -- MAIN MODELS (H2 / H3 / H4)
#
# H2: does exposure predict tree density / tree count?
# H3: does that relationship change with years since registration?
# H4: does the spatial concentration of exposure (10 thresholds) add to it?
#
# Farm area (Area_Ha_sc) is included in the TREE-COUNT models only.
# Density = Trees / Area_Ha, so adding farm area to the density model would
# use part of the outcome's own denominator as a predictor.
################################################################################

fit_density <- function(dat, nearfar = NULL) {
  lmer(build_formula("Density_winsor99", nearfar = nearfar), data = dat, REML = TRUE)
}
fit_trees <- function(dat, nearfar = NULL, area = TRUE, family = nbinom2) {
  glmmTMB(build_formula("Trees", area = area, nearfar = nearfar), family = family, data = dat)
}

fit_base_models <- function(dat) list(density = fit_density(dat), trees_nb = fit_trees(dat))
fit_h4_models   <- function(dat, nearfar_var) {
  list(density = fit_density(dat, nearfar_var), trees_nb = fit_trees(dat, nearfar_var))
}

message("Fitting base models (H2/H3)...")
base_B <- fit_base_models(TistDat_B)
base_S <- fit_base_models(TistDat_S)

message("Fitting H4 models across 10 thresholds...")
h4_B <- setNames(lapply(scaled_nf_cols, function(v) fit_h4_models(TistDat_B, v)), threshold_labels)
h4_S <- setNames(lapply(scaled_nf_cols, function(v) fit_h4_models(TistDat_S, v)), threshold_labels)

## ---- Convergence ---------------------------------------------------------------
convergence_row <- function(mod_list, label) {
  lmer_msg <- mod_list$density@optinfo$conv$lme4$messages
  tibble(
    model = label,
    density_status  = if (is.null(lmer_msg)) "OK" else paste(lmer_msg, collapse = "; "),
    density_singular = isSingular(mod_list$density),
    trees_status    = if (mod_list$trees_nb$fit$convergence == 0) "OK" else
      paste("code", mod_list$trees_nb$fit$convergence),
    trees_pdHess    = mod_list$trees_nb$sdr$pdHess
  )
}
store_result("3.1 Convergence -- main and H4 models", bind_rows(
  convergence_row(base_B, "Base -- Bushenyi"),
  convergence_row(base_S, "Base -- Soroti"),
  imap_dfr(h4_B, ~ convergence_row(.x, paste("H4", .y, "-- Bushenyi"))),
  imap_dfr(h4_S, ~ convergence_row(.x, paste("H4", .y, "-- Soroti")))
))

## ---- Main-model outputs (Table S8, Table 3.3) ---------------------------------
store_result("3.2 Main models: fixed effects (Table S8)", bind_rows(
  tidy_fixed(base_B$density,  "Bushenyi density (LMM)"),
  tidy_fixed(base_B$trees_nb, "Bushenyi tree count (NB GLMM, with farm area)"),
  tidy_fixed(base_S$density,  "Soroti density (LMM)"),
  tidy_fixed(base_S$trees_nb, "Soroti tree count (NB GLMM, with farm area)")
), csv = "main_models_fixed_effects.csv")

model_fit_row <- function(m, label) {
  tibble(
    model = label,
    n_obs = nobs(m),
    AIC   = round(AIC(m), 1),
    REML_criterion = if (inherits(m, "glmmTMB")) NA_real_ else round(REMLcrit(m), 1),
    NB_theta       = if (inherits(m, "glmmTMB")) round(sigma(m), 3) else NA_real_
  )
}
store_result("3.3 Main models: fit statistics (Table S8)", bind_rows(
  model_fit_row(base_B$density,  "Bushenyi density"),
  model_fit_row(base_B$trees_nb, "Bushenyi tree count"),
  model_fit_row(base_S$density,  "Soroti density"),
  model_fit_row(base_S$trees_nb, "Soroti tree count")
), csv = "main_models_fit.csv")

store_result("3.4 Main models: random-effect variances (Table S8)", c(
  "--- Bushenyi density ---",   capture.output(print(VarCorr(base_B$density), comp = c("Variance", "Std.Dev."))),
  "--- Bushenyi tree count ---", capture.output(print(VarCorr(base_B$trees_nb))),
  "--- Soroti density ---",     capture.output(print(VarCorr(base_S$density), comp = c("Variance", "Std.Dev."))),
  "--- Soroti tree count ---",   capture.output(print(VarCorr(base_S$trees_nb)))
))

## Tree-count models with vs without farm area (Section 3.9.4.1:
## "controlling for farm area weakened ... 31% to 20%").
trees_no_area_B <- fit_trees(TistDat_B, area = FALSE)
trees_no_area_S <- fit_trees(TistDat_S, area = FALSE)
store_result("3.5 Tree-count models with vs without farm area", bind_rows(
  tidy_fixed(trees_no_area_B, "Bushenyi, WITHOUT farm area"),
  tidy_fixed(base_B$trees_nb, "Bushenyi, WITH farm area"),
  tidy_fixed(trees_no_area_S, "Soroti, WITHOUT farm area"),
  tidy_fixed(base_S$trees_nb, "Soroti, WITH farm area")
) %>% filter(term != "(Intercept)"), csv = "tree_count_farm_area_comparison.csv")

## ---- AIC: does the H4 spatial-concentration term improve the density model? ---
## Refitted with ML (REML = FALSE) so AIC is comparable. Positive delta = H4 better.
compare_aic_h4 <- function(dat, site_label) {
  aic_base <- AIC(lmer(build_formula("Density_winsor99"), data = dat, REML = FALSE))
  aic_h4 <- sapply(scaled_nf_cols, function(v)
    AIC(lmer(build_formula("Density_winsor99", nearfar = v), data = dat, REML = FALSE)))
  tibble(site = site_label, threshold = threshold_labels,
         AIC_base = round(aic_base, 2), AIC_h4 = round(aic_h4, 2),
         delta_AIC = round(aic_base - aic_h4, 2))
}
store_result("3.6 AIC comparison, density: base vs H4 (positive delta = H4 better)",
             bind_rows(compare_aic_h4(TistDat_B, "Bushenyi"), compare_aic_h4(TistDat_S, "Soroti")))


################################################################################
# PART 4 -- DIAGNOSTICS AND FIGURES FOR THE MAIN MODELS
################################################################################

## ---- Collinearity between exposure and spatial concentration ------------------
## Correlation of the two fixed-effect estimates. |r| > 0.5 = the
## spatial-concentration coefficient is not reported as confirmed.
## Works for H4 lists ($density / $trees_nb, via `outcome`) and for the
## single-model H4b lists in Part 5 (outcome = NULL).
extract_collinearity <- function(mod_list, site_label, outcome = NULL) {
  map2_dfr(names(mod_list), scaled_nf_cols, function(th, nf_var) {
    m <- if (is.null(outcome)) mod_list[[th]] else mod_list[[th]][[outcome]]
    corr_mat <- cov2cor(get_vcov_matrix(m))
    tibble(threshold = th, site = site_label,
           exposure_nearfar_corr = round(corr_mat["Exposure_sc", nf_var], 3))
  })
}

## Coefficients of interest at every threshold. Near-Far term names change
## with the threshold, so they are relabelled to a common name.
extract_threshold_terms <- function(mod_list, site_label, outcome = NULL,
                                    base_terms = character(0), include_nearfar = TRUE) {
  map2_dfr(names(mod_list), scaled_nf_cols, function(th, nf_var) {
    m <- if (is.null(outcome)) mod_list[[th]] else mod_list[[th]][[outcome]]
    int_term <- paste0("Exposure_sc:", nf_var)
    keep <- c(base_terms, if (include_nearfar) c(nf_var, int_term))
    tidy_fixed(m, site_label) %>%
      filter(term %in% keep) %>%
      mutate(term = case_when(term == nf_var   ~ "NearFar_resid_sc",
                              term == int_term ~ "Exposure_sc:NearFar_resid_sc",
                              TRUE ~ term),
             threshold = th, site = site_label) %>%
      select(site, threshold, term, estimate, std.error, statistic, p.value)
  })
}

h4_base_terms <- c("Exposure_sc", "Exposure_sc:Years_since_reg_sc")

h4_collinearity_density <- bind_rows(extract_collinearity(h4_B, "Bushenyi", "density"),
                                     extract_collinearity(h4_S, "Soroti",   "density"))
h4_collinearity_trees   <- bind_rows(extract_collinearity(h4_B, "Bushenyi", "trees_nb"),
                                     extract_collinearity(h4_S, "Soroti",   "trees_nb"))

store_result("4.1 H4 own-outcome coefficients + collinearity -- density (Figure 3.12, 3.9.4.2-3)",
             bind_rows(extract_threshold_terms(h4_B, "Bushenyi", "density", h4_base_terms),
                       extract_threshold_terms(h4_S, "Soroti",   "density", h4_base_terms)) %>%
               left_join(h4_collinearity_density, by = c("threshold", "site")),
             csv = "h4_own_outcome_density.csv")
store_result("4.2 H4 own-outcome coefficients + collinearity -- tree count (Figure 3.12, 3.9.4.3)",
             bind_rows(extract_threshold_terms(h4_B, "Bushenyi", "trees_nb", h4_base_terms),
                       extract_threshold_terms(h4_S, "Soroti",   "trees_nb", h4_base_terms)) %>%
               left_join(h4_collinearity_trees, by = c("threshold", "site")),
             csv = "h4_own_outcome_trees.csv")

## ---- Predictions for Figures 3.10 and 3.12 -----------------------------------
generate_predictions <- function(model, predictor, site, outcome, effect_name,
                                 threshold = NA_character_) {
  ## No bias correction: ggeffects cannot compute the random-effect variance
  ## for these models (it warns that corrected results are "not reliable"),
  ## so predictions are for a typical group (random effects at zero).
  pred <- ggpredict(model, terms = predictor) %>%
    as.data.frame() %>%
    mutate(Site = site, Outcome = outcome, Effect = effect_name, Threshold = threshold)
  if ("group" %in% names(pred)) rename(pred, Duration_level = group) else mutate(pred, Duration_level = NA_character_)
}

## Back-transformation constants (pooled dataset, matching the scaling).
years_mu       <- mean(TistDat_H$Years_since_reg, na.rm = TRUE)
years_sigma    <- sd(TistDat_H$Years_since_reg,   na.rm = TRUE)
exposure_mu    <- mean(TistDat_H$Exposure, na.rm = TRUE)
exposure_sigma <- sd(TistDat_H$Exposure,   na.rm = TRUE)
nearfar1_mu    <- mean(TistDat_H$NearFar_resid_1, na.rm = TRUE)
nearfar1_sigma <- sd(TistDat_H$NearFar_resid_1,   na.rm = TRUE)

predict_four <- function(models, predictor, effect_name, threshold = NA_character_) {
  bind_rows(
    generate_predictions(models$B$density,  predictor, "Bushenyi", "Density", effect_name, threshold),
    generate_predictions(models$S$density,  predictor, "Soroti",   "Density", effect_name, threshold),
    generate_predictions(models$B$trees_nb, predictor, "Bushenyi", "Trees",   effect_name, threshold),
    generate_predictions(models$S$trees_nb, predictor, "Soroti",   "Trees",   effect_name, threshold)
  )
}
base_models <- list(B = base_B, S = base_S)

fig1_predictions <- predict_four(base_models, "Exposure_sc", "Exposure") %>%
  mutate(x_raw = x * exposure_sigma + exposure_mu)
duration_main_predictions <- predict_four(base_models, "Years_since_reg_sc", "Programme duration") %>%
  mutate(x_raw = x * years_sigma + years_mu)
duration_interaction_predictions <- predict_four(
  base_models, c("Exposure_sc", "Years_since_reg_sc [meansd]"), "Exposure × Duration"
) %>% mutate(x_raw = x * exposure_sigma + exposure_mu)
## Spatial concentration at 500 m only -- the only threshold not clearly collinear.
nearfar_predictions <- predict_four(list(B = h4_B[["500m"]], S = h4_S[["500m"]]),
                                    "NearFar_resid_1_sc", "Exposure concentration", "500m") %>%
  mutate(x_raw = x * nearfar1_sigma + nearfar1_mu)

plot_effect <- function(prediction_data, xlab) {
  pred_df <- prediction_data %>%
    mutate(across(c(x_raw, predicted, conf.low, conf.high), as.numeric),
           Outcome = factor(Outcome, levels = c("Density", "Trees"))) %>%
    arrange(Outcome, Site, x_raw)
  ggplot(pred_df, aes(x = x_raw, y = predicted, colour = Site, group = Site)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Site), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~Outcome, scales = "free_y") +
    labs(x = xlab, y = "Predicted value", colour = "Site", fill = "Site") +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
}

plot_duration_interaction <- function(prediction_data) {
  pred_df <- prediction_data %>%
    mutate(x_raw = as.numeric(x_raw),
           Duration_level = as.numeric(as.character(Duration_level))) %>%
    drop_na(Duration_level) %>%
    group_by(Site, Outcome) %>%
    mutate(Duration_level = factor(case_when(
      Duration_level == min(Duration_level) ~ "Short duration (-1 SD)",
      Duration_level == max(Duration_level) ~ "Long duration (+1 SD)",
      TRUE ~ "Average duration"),
      levels = c("Short duration (-1 SD)", "Average duration", "Long duration (+1 SD)"))) %>%
    ungroup() %>%
    arrange(Outcome, Site, Duration_level, x_raw)
  ggplot(pred_df, aes(x = x_raw, y = predicted, colour = Site, linetype = Duration_level,
                      group = interaction(Site, Duration_level))) +
    geom_line(linewidth = 1) +
    facet_wrap(~Outcome, scales = "free_y") +
    labs(x = "Exposure", y = "Predicted value", colour = "Site", linetype = "Programme duration") +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
}

fig1 <- plot_effect(fig1_predictions, "Exposure")
fig2 <- plot_effect(duration_main_predictions, "Years since registration") /
  plot_duration_interaction(duration_interaction_predictions)
fig3 <- plot_effect(nearfar_predictions, "Exposure concentration")

ggsave(file.path(fig_dir, "Figure1_exposure_effect.png"),     fig1, width = 10, height = 6,  dpi = 300)
ggsave(file.path(fig_dir, "Figure2_programme_duration.png"),  fig2, width = 10, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "Figure3_spatial_scale.png"),       fig3, width = 12, height = 7,  dpi = 300)

## ---- Bushenyi density: singular fit check --------------------------------------
## The full structure gives zero variance for district and sub-county.
## Refit with those levels dropped to confirm the exposure estimate is unaffected.
base_density_simplified_B <- lmer(build_formula("Density_winsor99", simplified = TRUE),
                                  data = TistDat_B, REML = TRUE)
store_result("4.3 Bushenyi density: full vs simplified random effects", bind_rows(
  tidy_fixed(base_B$density,            "Full (Admin_Districts/Subcounty/Village_ID)"),
  tidy_fixed(base_density_simplified_B, "Simplified (Village_ID only)")
) %>%
  filter(term == "Exposure_sc") %>%
  mutate(singular = c(isSingular(base_B$density), isSingular(base_density_simplified_B)),
         AIC = round(c(AIC(base_B$density), AIC(base_density_simplified_B)), 2)))

## ---- Farm-size check: are farms smaller in high-exposure neighbourhoods? -------
lm_coef_table <- function(m, label) {
  cf <- coef(summary(m))
  tibble(model = label, term = rownames(cf), estimate = cf[, 1], std.error = cf[, 2],
         statistic = cf[, 3], p.value = cf[, 4], r.squared = summary(m)$r.squared)
}
store_result("4.4 Farm size vs exposure (lm: Area_Ha_sc ~ Exposure_sc)", bind_rows(
  lm_coef_table(lm(Area_Ha_sc ~ Exposure_sc, data = TistDat_B), "Bushenyi"),
  lm_coef_table(lm(Area_Ha_sc ~ Exposure_sc, data = TistDat_S), "Soroti")
))

## ---- Soroti: group composition and cell sizes by tenure (Figure S4) -----------
soroti_group_composition <- TistDat_S %>%
  group_by(Group_ID) %>%
  summarise(n_members = n(), Duration_Tercile = first(Duration_Tercile), .groups = "drop") %>%
  mutate(group_size_cat = ifelse(n_members == 1, "Single-member", "Multi-member"))

group_composition_table <- soroti_group_composition %>%
  count(Duration_Tercile, group_size_cat) %>%
  group_by(Duration_Tercile) %>%
  mutate(pct_within_tercile = round(100 * n / sum(n), 1)) %>%
  ungroup()
store_result("4.5 Soroti group size by tenure tercile", group_composition_table)

group_size_chisq_tab <- table(soroti_group_composition$Duration_Tercile,
                              soroti_group_composition$group_size_cat)
chisq_result <- suppressWarnings(chisq.test(group_size_chisq_tab))
store_result("4.6 Chi-square: group size independent of tenure tercile (Soroti)",
             htest_row(chisq_result, "Soroti group size x tenure tercile"))

store_result("4.7 Soroti sample sizes: exposure tercile x tenure tercile", TistDat_S %>%
               mutate(Exposure_Tercile = factor(ntile(Exposure_sc, 3), labels = c("Low", "Medium", "High"))) %>%
               count(Exposure_Tercile, Duration_Tercile) %>%
               pivot_wider(names_from = Duration_Tercile, values_from = n, values_fill = 0))

tercile_labels <- TistDat_S %>%
  group_by(Duration_Tercile) %>%
  summarise(min_years = min(Years_since_reg, na.rm = TRUE),
            max_years = max(Years_since_reg, na.rm = TRUE), .groups = "drop") %>%
  mutate(year_label = ifelse(min_years == max_years, paste0(min_years, " yrs"),
                             paste0(min_years, "\u2013", max_years, " yrs")),
         axis_label = paste0(Duration_Tercile, " tenure\n(", year_label, ")")) %>%
  select(Duration_Tercile, axis_label)

plot_data <- group_composition_table %>%
  left_join(tercile_labels, by = "Duration_Tercile") %>%
  mutate(group_size_cat = factor(group_size_cat, levels = c("Multi-member", "Single-member")),
         axis_label = factor(axis_label, levels = tercile_labels$axis_label))

COL_MULTI  <- "#4C72B0"
COL_SINGLE <- "#DD8452"
fig_s4 <- ggplot(plot_data, aes(x = axis_label, y = pct_within_tercile, fill = group_size_cat)) +
  geom_col(width = 0.55, colour = NA, position = position_stack(reverse = TRUE)) +
  geom_text(data = filter(plot_data, group_size_cat == "Multi-member"),
            aes(label = paste0("n=", n)),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            colour = "white", fontface = "bold", size = 3.6) +
  geom_text(data = filter(plot_data, group_size_cat == "Single-member"),
            aes(y = 102, label = paste0(pct_within_tercile, "%\n(n=", n, ")")),
            colour = COL_SINGLE, fontface = "bold", size = 3.1, vjust = 0) +
  scale_fill_manual(values = c("Multi-member" = COL_MULTI, "Single-member" = COL_SINGLE)) +
  scale_y_continuous(limits = c(0, 113), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Share of Soroti groups (%)", fill = NULL,
       title = "Group composition by tenure tercile",
       subtitle = paste0("\u03c7\u00b2(", chisq_result$parameter, ") = ", round(chisq_result$statistic, 2),
                         ", p = ", format(round(chisq_result$p.value, 3), nsmall = 3),
                         if (chisq_result$p.value >= 0.05) " \u2014 not significant" else "")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic", colour = "grey30"),
        axis.text.x   = element_text(size = 10, lineheight = 0.9),
        axis.title.y  = element_text(size = 10.5))
ggsave(file.path(supp_fig_dir, "Figure_S4_soroti_tenure_composition.png"), fig_s4,
       width = 6.2, height = 5.6, dpi = 300)

## ---- Soroti tenure effect on tree count in plain terms (Figure 3.10 caption) ---
tenure_coef_S <- summary(base_S$trees_nb)$coefficients$cond["Years_since_reg_sc", "Estimate"]
store_result("4.8 Soroti tenure effect on tree count", tibble(
  coefficient = tenure_coef_S,
  pct_change_per_SD = 100 * (exp(tenure_coef_S) - 1)
))

tenure_tercile_means_S <- TistDat_S %>%
  group_by(Duration_Tercile) %>%
  summarise(mean_Years_since_reg_sc = mean(Years_since_reg_sc, na.rm = TRUE), .groups = "drop")

## Typical group (random effects at zero), other predictors at sample means;
## not bias-corrected (see generate_predictions()).
store_result("4.9 Soroti model-predicted tree count by tenure tercile (typical group)",
             as.data.frame(ggpredict(
               base_S$trees_nb,
               terms = paste0("Years_since_reg_sc [",
                              paste(round(tenure_tercile_means_S$mean_Years_since_reg_sc, 3), collapse = ","), "]")
             )) %>%
               mutate(Duration_Tercile = tenure_tercile_means_S$Duration_Tercile) %>%
               select(Duration_Tercile, x, predicted, conf.low, conf.high))

store_result("4.10 Soroti raw (unadjusted) tree count by tenure tercile", TistDat_S %>%
               group_by(Duration_Tercile) %>%
               summarise(n = n(), mean_Trees = round(mean(Trees, na.rm = TRUE), 1),
                         median_Trees = round(median(Trees, na.rm = TRUE), 1), .groups = "drop"))


################################################################################
# PART 5 -- NEIGHBOUR-DISSIMILARITY MODELS (H2b / H3b / H4b)
#
# Outcome = absolute difference between a farm and its neighbours. A
# NEGATIVE coefficient means neighbours become MORE alike as the predictor
# increases. Farm area is included for the tree-count outcome only (the
# density outcome is itself built from Trees / Area_Ha).
################################################################################

fit_dissim <- function(dat, dv_col, nearfar = NULL, area = FALSE, use_simplified = FALSE) {
  lmer(build_formula(dv_col, area = area, nearfar = nearfar, simplified = use_simplified),
       data = dat, REML = TRUE, na.action = na.exclude)
}

## Log-transformed tree-count dissimilarity columns.
for (k in 1:10) {
  TistDat_B[[paste0("log_AbsDissimilarity_Trees_", k)]] <- log1p(TistDat_B[[paste0("AbsDissimilarity_Trees_", k)]])
  TistDat_S[[paste0("log_AbsDissimilarity_Trees_", k)]] <- log1p(TistDat_S[[paste0("AbsDissimilarity_Trees_", k)]])
}
dv_cols_density <- paste0("AbsDissimilarity_Density_winsor99_", 1:10)
dv_cols_trees   <- paste0("log_AbsDissimilarity_Trees_", 1:10)

## Fits base_b (no Near-Far term) and h4b (with it) at all 10 thresholds.
fit_dissim_set <- function(dat, site, outcome, dv_cols, area) {
  simp <- map_lgl(threshold_labels, ~ use_simplified_re(site, outcome, .x))
  list(
    base = setNames(pmap(list(dv_cols, simp),
                         ~ fit_dissim(dat, ..1, area = area, use_simplified = ..2)), threshold_labels),
    h4b  = setNames(pmap(list(dv_cols, scaled_nf_cols, simp),
                         ~ fit_dissim(dat, ..1, nearfar = ..2, area = area, use_simplified = ..3)), threshold_labels)
  )
}

message("Fitting dissimilarity models (40 per outcome)...")
dis_density_B <- fit_dissim_set(TistDat_B, "Bushenyi", "Density", dv_cols_density, area = FALSE)
dis_density_S <- fit_dissim_set(TistDat_S, "Soroti",   "Density", dv_cols_density, area = FALSE)
dis_trees_B   <- fit_dissim_set(TistDat_B, "Bushenyi", "Trees",   dv_cols_trees,   area = TRUE)
dis_trees_S   <- fit_dissim_set(TistDat_S, "Soroti",   "Trees",   dv_cols_trees,   area = TRUE)

dissim_convergence <- function(mod_list, label) {
  imap_dfr(mod_list, function(m, th) {
    msg <- m@optinfo$conv$lme4$messages
    tibble(model = label, threshold = th,
           status = if (is.null(msg)) "OK" else paste(msg, collapse = "; "),
           singular = isSingular(m))
  })
}
store_result("5.1 Convergence -- dissimilarity models", bind_rows(
  dissim_convergence(dis_density_B$base, "base_b Density -- Bushenyi"),
  dissim_convergence(dis_density_S$base, "base_b Density -- Soroti"),
  dissim_convergence(dis_density_B$h4b,  "h4b Density -- Bushenyi"),
  dissim_convergence(dis_density_S$h4b,  "h4b Density -- Soroti"),
  dissim_convergence(dis_trees_B$base,   "base_b Trees -- Bushenyi"),
  dissim_convergence(dis_trees_S$base,   "base_b Trees -- Soroti"),
  dissim_convergence(dis_trees_B$h4b,    "h4b Trees -- Bushenyi"),
  dissim_convergence(dis_trees_S$h4b,    "h4b Trees -- Soroti")
))

compare_aic_b <- function(set, site_label, outcome) {
  tibble(site = site_label, outcome = outcome, threshold = threshold_labels,
         AIC_base_b = map_dbl(set$base, AIC), AIC_h4b = map_dbl(set$h4b, AIC)) %>%
    mutate(delta_AIC = round(AIC_base_b - AIC_h4b, 2))
}
store_result("5.2 AIC comparison, dissimilarity: base_b vs h4b (positive delta = h4b better)", bind_rows(
  compare_aic_b(dis_density_B, "Bushenyi", "Density"),
  compare_aic_b(dis_density_S, "Soroti",   "Density"),
  compare_aic_b(dis_trees_B,   "Bushenyi", "Trees"),
  compare_aic_b(dis_trees_S,   "Soroti",   "Trees")
))

## H2b / H3b: exposure and exposure x tenure on dissimilarity (Figures 3.9, 3.11)
h2b_h3b_density_table <- bind_rows(
  extract_threshold_terms(dis_density_B$base, "Bushenyi", base_terms = h4_base_terms, include_nearfar = FALSE),
  extract_threshold_terms(dis_density_S$base, "Soroti",   base_terms = h4_base_terms, include_nearfar = FALSE))
h2b_h3b_trees_table <- bind_rows(
  extract_threshold_terms(dis_trees_B$base, "Bushenyi", base_terms = h4_base_terms, include_nearfar = FALSE),
  extract_threshold_terms(dis_trees_S$base, "Soroti",   base_terms = h4_base_terms, include_nearfar = FALSE))
store_result("5.3 H2b/H3b coefficients -- density dissimilarity", h2b_h3b_density_table,
             csv = "h2b_h3b_density_dissimilarity.csv")
store_result("5.4 H2b/H3b coefficients -- tree-count dissimilarity (log1p)", h2b_h3b_trees_table,
             csv = "h2b_h3b_trees_dissimilarity.csv")

## H4b: spatial concentration on dissimilarity, with collinearity (Figure 3.13)
h4b_collinearity_density <- bind_rows(extract_collinearity(dis_density_B$h4b, "Bushenyi"),
                                      extract_collinearity(dis_density_S$h4b, "Soroti"))
h4b_collinearity_trees   <- bind_rows(extract_collinearity(dis_trees_B$h4b, "Bushenyi"),
                                      extract_collinearity(dis_trees_S$h4b, "Soroti"))

h4b_density_table <- bind_rows(extract_threshold_terms(dis_density_B$h4b, "Bushenyi"),
                               extract_threshold_terms(dis_density_S$h4b, "Soroti"))
h4b_trees_table   <- bind_rows(extract_threshold_terms(dis_trees_B$h4b, "Bushenyi"),
                               extract_threshold_terms(dis_trees_S$h4b, "Soroti"))

store_result("5.5 H4b coefficients + collinearity -- density dissimilarity",
             left_join(h4b_density_table, h4b_collinearity_density, by = c("threshold", "site")),
             csv = "h4b_density_dissimilarity.csv")
store_result("5.6 H4b coefficients + collinearity -- tree-count dissimilarity (log1p)",
             left_join(h4b_trees_table, h4b_collinearity_trees, by = c("threshold", "site")),
             csv = "h4b_trees_dissimilarity.csv")

## Collinearity panel, own-outcome and dissimilarity models (Figure S3).
collinearity_panel_df <- bind_rows(
  h4_collinearity_density  %>% mutate(outcome = "Density", model = "own-outcome"),
  h4_collinearity_trees    %>% mutate(outcome = "Trees",   model = "own-outcome"),
  h4b_collinearity_density %>% mutate(outcome = "Density", model = "dissimilarity"),
  h4b_collinearity_trees   %>% mutate(outcome = "Trees",   model = "dissimilarity")
) %>%
  mutate(threshold = factor(threshold, levels = threshold_labels),
         outcome   = factor(outcome, levels = c("Density", "Trees")),
         model     = factor(model, levels = c("own-outcome", "dissimilarity")),
         flagged   = abs(exposure_nearfar_corr) > 0.5)
store_result("5.7 Collinearity summary (Figure S3): all models, all thresholds",
             collinearity_panel_df, csv = "collinearity_all_models.csv")

fig_collinearity <- ggplot(collinearity_panel_df,
                           aes(x = threshold, y = exposure_nearfar_corr, colour = site, group = site)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = -0.5, fill = "grey70", alpha = 0.15) +
  geom_hline(yintercept = -0.5, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_grid(model ~ outcome) +
  scale_colour_manual(values = c("Bushenyi" = "#2E6F95", "Soroti" = "#C6602D")) +
  scale_y_continuous(limits = c(-0.9, 0), breaks = seq(-0.9, 0, 0.2)) +
  labs(x = "Exposure buffer threshold", y = "corr(Exposure_sc, NearFar_resid_sc)", colour = "Site",
       caption = "Shaded region: |r| > 0.5 -- coefficient flagged, not reported as confirmed") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        plot.caption = element_text(size = 8, colour = "grey40", face = "italic"),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"))
ggsave(file.path(supp_fig_dir, "h4_h4b_collinearity_panel.png"), fig_collinearity,
       width = 9, height = 7.5, dpi = 300)


################################################################################
# PART 6 -- COEFFICIENT-BY-THRESHOLD FIGURES (Figures 3.9, 3.11, 3.13)
################################################################################

prep_coef_table <- function(tbl, outcome_label) {
  tbl %>% mutate(Outcome = outcome_label,
                 threshold = factor(threshold, levels = threshold_labels),
                 ci_low = estimate - 1.96 * std.error,
                 ci_high = estimate + 1.96 * std.error)
}

plot_coef_by_threshold <- function(coef_data, facet_term = FALSE) {
  p <- ggplot(coef_data, aes(x = threshold, y = estimate, colour = site)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.4)) +
    labs(x = "Distance threshold", y = "Coefficient estimate (95% CI)", colour = "Site") +
    theme_classic(base_size = 15) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 13), legend.title = element_text(size = 13),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12), axis.title = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 13))
  if (facet_term) p + facet_grid(Outcome ~ term, scales = "free_y") else p + facet_wrap(~Outcome, scales = "free_y")
}

h2b_h3b_combined <- bind_rows(prep_coef_table(h2b_h3b_density_table, "Density"),
                              prep_coef_table(h2b_h3b_trees_table,   "Trees"))

fig_h2b_exposure    <- plot_coef_by_threshold(filter(h2b_h3b_combined, term == "Exposure_sc"))
fig_h3b_interaction <- plot_coef_by_threshold(filter(h2b_h3b_combined, term == "Exposure_sc:Years_since_reg_sc"))
fig_h4b <- bind_rows(prep_coef_table(h4b_density_table, "Density"),
                     prep_coef_table(h4b_trees_table,   "Trees")) %>%
  mutate(term = factor(term, levels = c("NearFar_resid_sc", "Exposure_sc:NearFar_resid_sc"),
                       labels = c("Near-Far residual (main effect)", "Exposure x Near-Far residual"))) %>%
  plot_coef_by_threshold(facet_term = TRUE)
fig_dissimilarity_combined <- (fig_h2b_exposure / fig_h3b_interaction) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figure_H2b_exposure_dissimilarity.png"), fig_h2b_exposure,    width = 11, height = 6,  dpi = 300)
ggsave(file.path(fig_dir, "figure_H3b_duration_dissimilarity.png"), fig_h3b_interaction, width = 11, height = 6,  dpi = 300)
ggsave(file.path(fig_dir, "figure_H4b_nearfar_dissimilarity.png"),  fig_h4b,             width = 13, height = 7,  dpi = 300)
ggsave(file.path(fig_dir, "figure_dissimilarity_combined.png"),     fig_dissimilarity_combined, width = 11, height = 11, dpi = 300)


################################################################################
# PART 7 -- NORM-CONVERGENCE ASYMMETRY TEST (Section 3.9.4.1)
#
# Signed dissimilarity (own - neighbour mean), fitted separately for farmers
# above and below their neighbourhood mean.
#   Above-mean: convergence predicts a NEGATIVE exposure coefficient.
#   Below-mean: convergence predicts a POSITIVE exposure coefficient.
################################################################################

test_convergence_asymmetry <- function(dat, k, site_label, use_simplified = FALSE) {
  signed_col <- paste0("Dissimilarity_Trees_", k)
  frm <- build_formula(signed_col, area = TRUE, simplified = use_simplified)
  
  map_dfr(c(TRUE, FALSE), function(above) {
    sub <- dat %>% filter((.data[[signed_col]] > 0) == above)
    group_label <- if (above) "Above neighbour mean" else "Below neighbour mean"
    m <- if (nrow(sub) >= 30) tryCatch(lmer(frm, data = sub, REML = TRUE, na.action = na.exclude),
                                       error = function(e) NULL) else NULL
    if (is.null(m)) {
      return(tibble(site = site_label, threshold = threshold_labels[k], group = group_label,
                    n = nrow(sub), exposure_est = NA_real_, exposure_se = NA_real_,
                    exposure_t = NA_real_, exposure_p = NA_real_))
    }
    cf <- tidy_fixed(m, "") %>% filter(term == "Exposure_sc")
    tibble(site = site_label, threshold = threshold_labels[k], group = group_label, n = nrow(sub),
           exposure_est = cf$estimate, exposure_se = cf$std.error,
           exposure_t = cf$statistic, exposure_p = cf$p.value)
  })
}

convergence_thresholds <- c(2, 4, 6, 8)   # 1,000 / 2,000 / 3,000 / 4,000 m
store_result("7.1 Norm-convergence asymmetry: exposure on signed tree-count dissimilarity", bind_rows(
  map_dfr(convergence_thresholds, ~ test_convergence_asymmetry(
    TistDat_B, .x, "Bushenyi", use_simplified_re("Bushenyi", "Trees", threshold_labels[.x]))),
  map_dfr(convergence_thresholds, ~ test_convergence_asymmetry(
    TistDat_S, .x, "Soroti", use_simplified_re("Soroti", "Trees", threshold_labels[.x])))
), csv = "convergence_asymmetry.csv")


################################################################################
# PART 8 -- RESIDUAL DIAGNOSTICS FOR THE MAIN MODELS (Section 3.9.3, Table S4)
#
#   Moran's I: village-mean residuals, 4 nearest neighbours.
#   DHARMa:    simulation-based residuals at farmer level; spatial
#              autocorrelation, dispersion, zero-inflation.
#   Zero-truncated NB: every farmer planted at least one tree, so tree
#              counts contain no zeros; the tree-count models are refitted
#              with a zero-truncated negative binomial as a robustness check.
################################################################################

## ---- Moran's I on village-aggregated residuals --------------------------------
dat_B <- TistDat_B %>% mutate(res_density = resid(base_B$density),
                              res_trees   = residuals(base_B$trees_nb, type = "pearson"))
dat_S <- TistDat_S %>% mutate(res_density = resid(base_S$density),
                              res_trees   = residuals(base_S$trees_nb, type = "pearson"))

village_residuals <- function(dat) {
  dat %>% group_by(Village_ID) %>%
    summarise(res_density = mean(res_density, na.rm = TRUE), res_trees = mean(res_trees, na.rm = TRUE),
              lon = mean(longitude, na.rm = TRUE), lat = mean(latitude, na.rm = TRUE), .groups = "drop")
}
village_res_B <- village_residuals(dat_B)
village_res_S <- village_residuals(dat_S)

build_spatial_weights <- function(village_res, utm_crs = 32736, k = 4) {
  coords <- st_as_sf(village_res, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(utm_crs) %>% st_coordinates()
  nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")
}
lw_B <- build_spatial_weights(village_res_B)
lw_S <- build_spatial_weights(village_res_S)

store_result("8.1 Moran's I on village-mean residuals", bind_rows(
  htest_row(moran.test(village_res_B$res_density, lw_B), "Bushenyi density"),
  htest_row(moran.test(village_res_B$res_trees,   lw_B), "Bushenyi tree count (Pearson)"),
  htest_row(moran.test(village_res_S$res_density, lw_S), "Soroti density"),
  htest_row(moran.test(village_res_S$res_trees,   lw_S), "Soroti tree count (Pearson)")
), csv = "morans_i.csv")

## ---- DHARMa ----------------------------------------------------------------------
## Duplicated coordinates (several records at one grove location) make
## testSpatialAutocorrelation() fail, so residuals are aggregated by location
## first when duplicates exist.
test_spatial_dharma <- function(sim, dat) {
  loc <- factor(paste(dat$longitude, dat$latitude))
  if (anyDuplicated(loc) > 0) {
    sim_use <- recalculateResiduals(sim, group = loc)
    coords  <- tibble(loc = loc, x = dat$longitude, y = dat$latitude) %>%
      group_by(loc) %>% summarise(x = first(x), y = first(y), .groups = "drop")
  } else {
    sim_use <- sim
    coords  <- tibble(x = dat$longitude, y = dat$latitude)
  }
  testSpatialAutocorrelation(sim_use, x = coords$x, y = coords$y, plot = FALSE)
}

message("Simulating DHARMa residuals (", n_dharma_sim, " simulations x 4 models)...")
sims <- list(
  "Bushenyi density"    = simulateResiduals(base_B$density,  n = n_dharma_sim),
  "Bushenyi tree count" = simulateResiduals(base_B$trees_nb, n = n_dharma_sim),
  "Soroti density"      = simulateResiduals(base_S$density,  n = n_dharma_sim),
  "Soroti tree count"   = simulateResiduals(base_S$trees_nb, n = n_dharma_sim)
)
sim_data <- list(dat_B, dat_B, dat_S, dat_S)

if (save_dharma_plots) {
  iwalk(sims, function(s, label) {
    png(file.path(supp_fig_dir, paste0("DHARMa_", gsub(" ", "_", label), ".png")),
        width = 2400, height = 1200, res = 200)
    plot(s, main = label)
    dev.off()
  })
}

store_result("8.2 DHARMa spatial autocorrelation (farmer-level residuals)",
             bind_rows(map2(sims, sim_data, test_spatial_dharma) %>%
                         imap(~ htest_row(.x, .y))),
             csv = "dharma_spatial.csv")

store_result("8.3 DHARMa dispersion and zero-inflation (tree-count models)", bind_rows(
  htest_row(testDispersion(sims[["Bushenyi tree count"]],    plot = FALSE), "Dispersion -- Bushenyi tree count"),
  htest_row(testDispersion(sims[["Soroti tree count"]],      plot = FALSE), "Dispersion -- Soroti tree count"),
  htest_row(testZeroInflation(sims[["Bushenyi tree count"]], plot = FALSE), "Zero-inflation -- Bushenyi tree count"),
  htest_row(testZeroInflation(sims[["Soroti tree count"]],   plot = FALSE), "Zero-inflation -- Soroti tree count")
), csv = "dharma_dispersion_zeroinflation.csv")

## ---- Zero-truncated negative binomial robustness check --------------------------
## Same specification as the main tree-count models (build_formula), with a
## zero-truncated family. Key terms should keep their direction and significance.
message("Fitting zero-truncated tree-count models...")
trunc_B <- fit_trees(TistDat_B, family = truncated_nbinom2)
trunc_S <- fit_trees(TistDat_S, family = truncated_nbinom2)

store_result("8.4 Zero-truncated NB robustness check (Table S4)", bind_rows(
  tidy_fixed(base_B$trees_nb, "Bushenyi, NB"),
  tidy_fixed(trunc_B,         "Bushenyi, truncated NB"),
  tidy_fixed(base_S$trees_nb, "Soroti, NB"),
  tidy_fixed(trunc_S,         "Soroti, truncated NB")
) %>% filter(term != "(Intercept)") %>% arrange(term, model),
csv = "zero_truncated_comparison.csv")

store_result("8.5 Zero-truncated NB convergence", tibble(
  model  = c("Bushenyi, truncated NB", "Soroti, truncated NB"),
  status = c(trunc_B$fit$convergence, trunc_S$fit$convergence),
  pdHess = c(trunc_B$sdr$pdHess, trunc_S$sdr$pdHess)
))


################################################################################
# PART 9 -- RESULTS SUMMARY
#
# Prints every stored result, in order, and writes the same output to
# <results_dir>/results_summary.txt. Tables are also saved as CSV above.
################################################################################

print_result <- function(name, obj) {
  cat("\n", strrep("=", 80), "\n", name, "\n", strrep("=", 80), "\n", sep = "")
  if (is.character(obj)) {
    cat(obj, sep = "\n")
  } else if (inherits(obj, "tbl_df")) {
    print(obj, n = Inf, width = Inf)
  } else {
    print(obj)
  }
}

summary_file <- file.path(results_dir, "results_summary.txt")
sink(summary_file, split = TRUE)
cat("TIST MODELLING PIPELINE -- RESULTS SUMMARY\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M"), "| R", as.character(getRversion()),
    "| lme4", as.character(packageVersion("lme4")),
    "| glmmTMB", as.character(packageVersion("glmmTMB")),
    "| DHARMa", as.character(packageVersion("DHARMa")), "\n")
cat("Settings: variogram_complete_cases =", variogram_complete_cases,
    "| n_dharma_sim =", n_dharma_sim, "| seed = 2026\n")
iwalk(results, ~ print_result(.y, .x))
cat("\nSignificance guide: lmer p-values use the normal approximation to t (|t| > 2 ~ p < .05).\n")
cat("Collinearity guide: |exposure_nearfar_corr| > 0.5 = not reported as confirmed.\n")
cat("Outputs written to:", normalizePath(output_dir), "\n")
sink()

message("Done. Full results summary: ", summary_file)

################################################################################
# END OF SCRIPT
################################################################################