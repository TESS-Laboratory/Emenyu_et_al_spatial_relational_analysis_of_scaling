# ============================================================
# TIST MIXED-MODEL / VARIOGRAM / SPATIAL-DIAGNOSTICS FUNCTIONS
#
# Continuation of pipeline_functions_v2.R. These sections did
# NOT exist in the earlier (v1) function library -- they are new,
# refactored directly from the tested flat script's mixed-model,
# similarity-model, variogram, and spatial-diagnostic code.
# ============================================================


# ============================================================
# SECTION F: PART 0/1 -- SETUP + MAIN MODELS (H2/H3/H4)
# ============================================================

get_scaled_nf_cols <- function() {
  paste0("NearFar_resid_", 1:10, "_sc")
}


get_threshold_labels <- function() {
  paste0(seq(500, 5000, by = 500), "m")
}


# ------------------------------------------------------------
# Add all scaled predictors used throughout the modelling
# sections, matching the tested script's Part 0 mutate().
# ------------------------------------------------------------

add_scaled_predictors <- function(dat) {
  
  dat |>
    dplyr::mutate(
      Dist_To_Forest_sc  = as.numeric(scale(Dist_To_Forest_m)),
      Area_Ha_sc         = as.numeric(scale(Area_Ha)),
      Trees_sc           = as.numeric(scale(Trees)),
      Exposure_sc        = as.numeric(scale(Exposure)),
      Years_since_reg_sc = as.numeric(scale(Years_since_reg)),
      dplyr::across(
        NearFar_resid_1:NearFar_resid_10,
        ~ as.numeric(scale(.x)),
        .names = "{.col}_sc"
      )
    )
}


# ------------------------------------------------------------
# Diagnostic: correlation of each scaled Near-Far column with
# Exposure_sc (should be ~0 by construction).
# ------------------------------------------------------------

check_nearfar_exposure_independence <- function(
    dat,
    scaled_nf_cols = get_scaled_nf_cols()
) {
  
  sapply(scaled_nf_cols, function(v) {
    round(
      cor(
        dat[[v]],
        dat$Exposure_sc,
        use = "complete.obs"
      ),
      8
    )
  })
}


# ------------------------------------------------------------
# Split data by study site
# ------------------------------------------------------------

split_by_site <- function(dat) {
  
  list(
    Bushenyi = dat |> dplyr::filter(Proj_Area == "Bushenyi"),
    Soroti   = dat |> dplyr::filter(Proj_Area == "Soroti")
  )
}


# ------------------------------------------------------------
# Splits a site's farms into Short/Medium/Long tenure terciles
# based on Years_since_reg, cut separately per site.
# ------------------------------------------------------------

build_duration_tercile <- function(dat) {
  
  dat |>
    dplyr::group_by(Proj_Area) |>
    dplyr::mutate(
      Duration_Tercile = dplyr::case_when(
        
        Years_since_reg <= quantile(
          Years_since_reg,
          1 / 3,
          na.rm = TRUE
        ) ~ "Short",
        
        Years_since_reg <= quantile(
          Years_since_reg,
          2 / 3,
          na.rm = TRUE
        ) ~ "Medium",
        
        TRUE ~ "Long"
      ),
      
      Duration_Tercile = factor(
        Duration_Tercile,
        levels = c(
          "Short",
          "Medium",
          "Long"
        )
      )
    ) |>
    dplyr::ungroup()
}


# ------------------------------------------------------------
# H2/H3 base models:
# density = lmer
# tree count = glmmTMB negative-binomial nbinom2
#
# NOTE:
# Area_Ha_sc is included in the trees formula only -- NOT in
# density, since Density_winsor99 = Trees / Area_Ha.
# Including Area_Ha_sc in the density model would therefore
# use part of the outcome's denominator to predict the outcome.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Base models
#
# Density outcome:
# Density_winsor99 is modelled with a simplified random-effects
# structure: Village_ID, Cluster_ID and Group_ID.
#
# Justification:
# The previous nested structure
# (1 | Admin_Districts/Subcounty/Village_ID)
# produced singular/near-zero variance components in the
# density models, particularly for Bushenyi. A previous
# diagnostic showed that removing the higher-level
# Admin_Districts/Subcounty components resolved the degenerate
# structure without materially changing the Exposure_sc effect.
#
# Village_ID is retained because it captures local spatial
# clustering, while Cluster_ID and Group_ID capture the
# programme's finer social/spatial organisation.
#
# Area_Ha_sc is deliberately NOT included in the density model
# because Density_winsor99 is derived from Trees / Area_Ha.
# Including farm area would therefore introduce a predictor
# that forms part of the outcome's denominator.
# ------------------------------------------------------------

fit_base_models <- function(dat) {
  
  list(
    
    density = lme4::lmer(
      Density_winsor99 ~
        Exposure_sc +
        Years_since_reg_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Village_ID) +
        (1 | Cluster_ID) +
        (1 | Group_ID),
      data = dat,
      REML = TRUE
    ),
    
    # Tree counts retain the original hierarchical structure
    # because Trees is a count outcome and farm size is a
    # legitimate exposure/scale variable for this outcome.
    trees_nb = glmmTMB::glmmTMB(
      Trees ~
        Exposure_sc +
        Years_since_reg_sc +
        Area_Ha_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Admin_Districts/Subcounty/Village_ID) +
        (1 | Cluster_ID/Group_ID),
      family = glmmTMB::nbinom2,
      data = dat
    )
  )
}
# ------------------------------------------------------------
# H4 extension:
# Base model + one Near-Far spatial-concentration term
# + its interaction with Exposure.
#
# Run once per threshold.
# ------------------------------------------------------------

fit_h4_models <- function(dat, nearfar_var) {
  density_formula <- as.formula(
    paste0(
      "Density_winsor99 ~ ",
      "Exposure_sc + ",
      "Years_since_reg_sc + ",
      "Exposure_sc:Years_since_reg_sc + ",
      "Dist_To_Forest_sc + ",
      nearfar_var,
      " + Exposure_sc:",
      nearfar_var,
      " + ",
      "(1 | Village_ID) + (1 | Cluster_ID) + (1 | Group_ID)"
    )
  )
  
  # trees_formula unchanged -- matches fit_base_models()'s trees_nb,
  # which keeps the original nested structure.
  trees_formula <- as.formula(
    paste0(
      "Trees ~ ",
      "Exposure_sc + ",
      "Years_since_reg_sc + ",
      "Area_Ha_sc + ",
      "Exposure_sc:Years_since_reg_sc + ",
      "Dist_To_Forest_sc + ",
      nearfar_var,
      " + Exposure_sc:",
      nearfar_var,
      " + ",
      "(1 | Admin_Districts/Subcounty/Village_ID) + ",
      "(1 | Cluster_ID/Group_ID)"
    )
  )
  
  list(
    density = lme4::lmer(
      density_formula,
      data = dat,
      REML = TRUE
    ),
    trees_nb = glmmTMB::glmmTMB(
      trees_formula,
      family = glmmTMB::nbinom2,
      data = dat
    )
  )
}

# ------------------------------------------------------------
# Fits fit_h4_models() across all 10 thresholds.
#
# Returns a named list with one element per threshold.
# Kept as ONE target per site rather than 10 separate targets,
# matching the tested script's setNames(lapply(...)) pattern.
# ------------------------------------------------------------

fit_h4_models_all_thresholds <- function(
    dat,
    scaled_nf_cols = get_scaled_nf_cols(),
    threshold_labels = get_threshold_labels()
) {
  
  setNames(
    lapply(
      scaled_nf_cols,
      function(v) {
        fit_h4_models(
          dat,
          nearfar_var = v
        )
      }
    ),
    threshold_labels
  )
}


# ------------------------------------------------------------
# Model convergence diagnostics
# ------------------------------------------------------------

check_convergence <- function(mod_list, label) {
  
  lmer_warn <- mod_list$density@optinfo$conv$lme4$messages
  
  nb_conv <- mod_list$trees_nb$fit$convergence
  
  tibble::tibble(
    label = label,
    
    density_convergence_ok = is.null(lmer_warn),
    
    density_singular =
      lme4::isSingular(
        mod_list$density,
        tol = 1e-4
      ),
    
    density_warning = ifelse(
      is.null(lmer_warn),
      NA_character_,
      paste(
        lmer_warn,
        collapse = "; "
      )
    ),
    
    trees_convergence_ok = nb_conv == 0,
    
    trees_convergence_code = nb_conv
  )
}


check_convergence_across_thresholds <- function(
    mod_list_by_threshold,
    label
) {
  
  purrr::map_dfr(
    names(mod_list_by_threshold),
    function(th) {
      
      check_convergence(
        mod_list_by_threshold[[th]],
        label
      ) |>
        dplyr::mutate(
          threshold = th,
          .after = label
        )
    }
  )
}


# ------------------------------------------------------------
# AIC comparison:
# Does adding the Near-Far term (H4) improve on the base
# density model at each of the 10 thresholds?
# ------------------------------------------------------------

compare_aic <- function(
    dat,
    scaled_nf_cols = get_scaled_nf_cols(),
    threshold_labels = get_threshold_labels(),
    site_label
) {
  
  base_ml <- lme4::lmer(
    Density_winsor99 ~
      Exposure_sc +
      Years_since_reg_sc +
      Exposure_sc:Years_since_reg_sc +
      Dist_To_Forest_sc +
      (1 | Village_ID) +
      (1 | Cluster_ID) +
      (1 | Group_ID),
    data = dat,
    REML = FALSE
  )
  
  aic_base <- AIC(base_ml)
  
  aic_h4 <- sapply(
    scaled_nf_cols,
    function(v) {
      
      f <- as.formula(
        paste0(
          "Density_winsor99 ~ ",
          "Exposure_sc + ",
          "Years_since_reg_sc + ",
          "Exposure_sc:Years_since_reg_sc + ",
          "Dist_To_Forest_sc + ",
          v,
          " + Exposure_sc:",
          v,
          " + ",
          "(1 | Village_ID) + (1 | Cluster_ID) + (1 | Group_ID)"
        )
      )
      
      AIC(
        lme4::lmer(
          f,
          data = dat,
          REML = FALSE
        )
      )
    }
  )
  
  data.frame(
    site = site_label,
    threshold = threshold_labels,
    AIC_base = round(aic_base, 2),
    AIC_h4 = round(aic_h4, 2),
    delta_AIC = round(
      aic_base - aic_h4,
      2
    )
  )
}

# ------------------------------------------------------------
# Generate model predictions using ggeffects
# ------------------------------------------------------------

generate_predictions <- function(
    model,
    predictor,
    site,
    outcome,
    effect_name,
    threshold = NA_character_
) {
  
  pred <- ggeffects::ggpredict(
    model,
    terms = predictor
  ) |>
    as.data.frame()
  
  pred <- pred |>
    dplyr::mutate(
      Site = site,
      Outcome = outcome,
      Effect = effect_name,
      Threshold = threshold
    )
  
  if ("group" %in% names(pred)) {
    
    pred <- pred |>
      dplyr::rename(
        Duration_level = group
      )
    
  } else {
    
    pred$Duration_level <- NA_character_
  }
  
  pred
}


# ------------------------------------------------------------
# Compute global scaling constants used to back-transform
# scaled predictors for plotting.
# ------------------------------------------------------------

compute_global_scaling_constants <- function(dat) {
  
  nearfar_mu_sigma <- lapply(
    1:10,
    function(k) {
      
      col <- dat[[paste0(
        "NearFar_resid_",
        k
      )]]
      
      c(
        mu = mean(
          col,
          na.rm = TRUE
        ),
        
        sigma = sd(
          col,
          na.rm = TRUE
        )
      )
    }
  )
  
  names(nearfar_mu_sigma) <- as.character(1:10)
  
  list(
    
    years_mu =
      mean(
        dat$Years_since_reg,
        na.rm = TRUE
      ),
    
    years_sigma =
      sd(
        dat$Years_since_reg,
        na.rm = TRUE
      ),
    
    exposure_mu =
      mean(
        dat$Exposure,
        na.rm = TRUE
      ),
    
    exposure_sigma =
      sd(
        dat$Exposure,
        na.rm = TRUE
      ),
    
    nearfar_mu_sigma =
      nearfar_mu_sigma
  )
}


# ------------------------------------------------------------
# Figure 1 prediction data
# ------------------------------------------------------------

build_fig1_predictions <- function(
    base_B,
    base_S,
    scaling
) {
  
  dplyr::bind_rows(
    
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
    
  ) |>
    dplyr::mutate(
      x_raw =
        x * scaling$exposure_sigma +
        scaling$exposure_mu
    )
}


# ------------------------------------------------------------
# Programme-duration main-effect predictions
# ------------------------------------------------------------

build_duration_main_predictions <- function(
    base_B,
    base_S,
    scaling
) {
  
  dplyr::bind_rows(
    
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
    
  ) |>
    dplyr::mutate(
      x_raw =
        x * scaling$years_sigma +
        scaling$years_mu
    )
}


# ------------------------------------------------------------
# Exposure × programme-duration interaction predictions
# ------------------------------------------------------------

build_duration_interaction_predictions <- function(
    base_B,
    base_S,
    scaling
) {
  
  term <- c(
    "Exposure_sc",
    "Years_since_reg_sc [meansd]"
  )
  
  dplyr::bind_rows(
    
    generate_predictions(
      base_B$density,
      term,
      "Bushenyi",
      "Density",
      "Exposure × Duration"
    ),
    
    generate_predictions(
      base_S$density,
      term,
      "Soroti",
      "Density",
      "Exposure × Duration"
    ),
    
    generate_predictions(
      base_B$trees_nb,
      term,
      "Bushenyi",
      "Trees",
      "Exposure × Duration"
    ),
    
    generate_predictions(
      base_S$trees_nb,
      term,
      "Soroti",
      "Trees",
      "Exposure × Duration"
    )
    
  ) |>
    dplyr::mutate(
      x_raw =
        x * scaling$exposure_sigma +
        scaling$exposure_mu
    )
}


# ------------------------------------------------------------
# Near-Far predictions
#
# Default plotting thresholds:
# 500m  = index 1
# 2000m = index 4
# 3500m = index 7
# ------------------------------------------------------------

build_nearfar_predictions <- function(
    h4_B,
    h4_S,
    scaling,
    thresholds = c(
      "500m" = 1,
      "2000m" = 4,
      "3500m" = 7
    )
) {
  
  nearfar_predictions <-
    purrr::imap_dfr(
      thresholds,
      function(index, threshold) {
        
        variable <- paste0(
          "NearFar_resid_",
          index,
          "_sc"
        )
        
        dplyr::bind_rows(
          
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
    ) |>
    dplyr::mutate(
      Threshold = factor(
        Threshold,
        levels = names(thresholds)
      )
    )
  
  nearfar_lookup <- tibble::tibble(
    
    Threshold = factor(
      names(thresholds),
      levels = names(thresholds)
    ),
    
    mu = vapply(
      thresholds,
      function(i) {
        scaling$nearfar_mu_sigma[[
          as.character(i)
        ]]["mu"]
      },
      numeric(1)
    ),
    
    sigma = vapply(
      thresholds,
      function(i) {
        scaling$nearfar_mu_sigma[[
          as.character(i)
        ]]["sigma"]
      },
      numeric(1)
    )
  )
  
  nearfar_predictions |>
    dplyr::left_join(
      nearfar_lookup,
      by = "Threshold"
    ) |>
    dplyr::mutate(
      x_raw = x * sigma + mu
    ) |>
    dplyr::select(
      -mu,
      -sigma
    )
}


# ------------------------------------------------------------
# Generic effect plot
# ------------------------------------------------------------

plot_effect <- function(
    prediction_data,
    xlab,
    facet_threshold = FALSE
) {
  
  pred_df <- prediction_data |>
    dplyr::mutate(
      x_raw = as.numeric(x_raw),
      predicted = as.numeric(predicted),
      conf.low = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      
      Outcome = factor(
        Outcome,
        levels = c(
          "Density",
          "Trees"
        )
      )
    ) |>
    dplyr::arrange(
      Outcome,
      Site,
      Threshold,
      x_raw
    ) |>
    dplyr::mutate(
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
  
  p <- ggplot2::ggplot(
    pred_df,
    ggplot2::aes(
      x = x_raw,
      y = predicted,
      colour = Site,
      group = Plot_group
    )
  ) +
    
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = conf.low,
        ymax = conf.high,
        fill = Site,
        group = Plot_group
      ),
      alpha = 0.12,
      colour = NA
    ) +
    
    ggplot2::geom_line(
      linewidth = 1
    )
  
  p <- if (facet_threshold) {
    
    p +
      ggplot2::facet_grid(
        Outcome ~ Threshold,
        scales = "free_y"
      )
    
  } else {
    
    p +
      ggplot2::facet_wrap(
        ~Outcome,
        scales = "free_y"
      )
  }
  
  p +
    ggplot2::labs(
      x = xlab,
      y = "Predicted value",
      colour = "Site",
      fill = "Site"
    ) +
    
    ggplot2::theme_classic(
      base_size = 13
    ) +
    
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    )
}


# ------------------------------------------------------------
# Plot Exposure × Duration interaction
# ------------------------------------------------------------

plot_duration_interaction <- function(
    prediction_data
) {
  
  pred_df <- prediction_data |>
    dplyr::filter(
      Effect == "Exposure × Duration"
    ) |>
    
    dplyr::mutate(
      x_raw = as.numeric(x_raw),
      Duration_level =
        as.numeric(
          as.character(Duration_level)
        )
    ) |>
    
    tidyr::drop_na(
      Duration_level
    ) |>
    
    dplyr::group_by(
      Site,
      Outcome
    ) |>
    
    dplyr::mutate(
      
      Duration_level =
        dplyr::case_when(
          
          Duration_level ==
            min(Duration_level) ~
            "Short duration (-1 SD)",
          
          Duration_level ==
            max(Duration_level) ~
            "Long duration (+1 SD)",
          
          TRUE ~
            "Average duration"
        ),
      
      Duration_level =
        factor(
          Duration_level,
          levels = c(
            "Short duration (-1 SD)",
            "Average duration",
            "Long duration (+1 SD)"
          )
        )
    ) |>
    
    dplyr::ungroup() |>
    
    dplyr::arrange(
      Outcome,
      Site,
      Duration_level,
      x_raw
    )
  
  ggplot2::ggplot(
    pred_df,
    ggplot2::aes(
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
    
    ggplot2::geom_line(
      linewidth = 1
    ) +
    
    ggplot2::facet_wrap(
      ~Outcome,
      scales = "free_y"
    ) +
    
    ggplot2::labs(
      x = "Exposure",
      y = "Predicted value",
      colour = "Site",
      linetype = "Programme duration"
    ) +
    
    ggplot2::theme_classic(
      base_size = 13
    ) +
    
    ggplot2::theme(
      legend.position = "bottom"
    )
}


# ------------------------------------------------------------
# Save ggplot
# ------------------------------------------------------------

save_ggplot <- function(
    plot,
    path,
    width,
    height,
    dpi = 300
) {
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggplot2::ggsave(
    path,
    plot,
    width = width,
    height = height,
    dpi = dpi
  )
  
  path
}


# ------------------------------------------------------------
# Bushenyi density diagnostic:
#
# Re-fits H2 only with the two zero-variance random-effect
# levels dropped, and compares Exposure_sc's estimate/SE/t
# between structures.
#
# IMPORTANT:
# The original supplied version had four values inside
# singular_fit for a two-row tibble. This corrected version
# keeps singularity and convergence as separate diagnostics.
# ------------------------------------------------------------

diagnose_bushenyi_density_simplified_re <- function(
    TistDat_B,
    base_B_density
) {
  
  base_density_simplified_B <- lme4::lmer(
    
    Density_winsor99 ~
      Exposure_sc +
      Years_since_reg_sc +
      Exposure_sc:Years_since_reg_sc +
      Dist_To_Forest_sc +
      (1 | Village_ID) +
      (1 | Cluster_ID/Group_ID),
    
    data = TistDat_B,
    REML = TRUE
  )
  
  diag_original <-
    coef(summary(base_B_density))[
      "Exposure_sc",
    ]
  
  diag_simplified <-
    coef(summary(base_density_simplified_B))[
      "Exposure_sc",
    ]
  
  tibble::tibble(
    
    structure = c(
      "Original (Admin_Districts/Subcounty/Village_ID)",
      "Simplified (Village_ID only)"
    ),
    
    singular_fit = c(
      lme4::isSingular(
        base_B_density,
        tol = 1e-4
      ),
      
      lme4::isSingular(
        base_density_simplified_B,
        tol = 1e-4
      )
    ),
    
    convergence_warning = c(
      
      !is.null(
        base_B_density@optinfo$conv$lme4$messages
      ),
      
      !is.null(
        base_density_simplified_B@optinfo$conv$lme4$messages
      )
    ),
    
    AIC = round(
      c(
        AIC(base_B_density),
        AIC(base_density_simplified_B)
      ),
      2
    ),
    
    Exposure_sc_estimate = round(
      c(
        diag_original["Estimate"],
        diag_simplified["Estimate"]
      ),
      2
    ),
    
    Exposure_sc_se = round(
      c(
        diag_original["Std. Error"],
        diag_simplified["Std. Error"]
      ),
      2
    ),
    
    Exposure_sc_t = round(
      c(
        diag_original["t value"],
        diag_simplified["t value"]
      ),
      3
    )
  )
}


# ------------------------------------------------------------
# Farm-size confound check:
# Does Exposure correlate with farm size?
# ------------------------------------------------------------

check_farm_size_confound <- function(
    dat,
    site_label
) {
  
  model <- lm(
    Area_Ha_sc ~ Exposure_sc,
    data = dat
  )
  
  list(
    site = site_label,
    model = model,
    summary = summary(model)
  )
}


# ------------------------------------------------------------
# Soroti check 1:
# Is single-member-group status independent of tenure tercile?
# ------------------------------------------------------------

check_soroti_group_composition <- function(
    TistDat_S
) {
  
  soroti_group_composition <-
    TistDat_S |>
    dplyr::group_by(
      Group_ID
    ) |>
    
    dplyr::summarise(
      
      n_members = dplyr::n(),
      
      Duration_Tercile =
        dplyr::first(
          Duration_Tercile
        ),
      
      .groups = "drop"
    ) |>
    
    dplyr::mutate(
      
      group_size_cat =
        ifelse(
          n_members == 1,
          "Single-member",
          "Multi-member"
        )
    )
  
  group_composition_table <-
    soroti_group_composition |>
    
    dplyr::count(
      Duration_Tercile,
      group_size_cat
    ) |>
    
    dplyr::group_by(
      Duration_Tercile
    ) |>
    
    dplyr::mutate(
      pct_within_tercile =
        round(
          100 * n / sum(n),
          1
        )
    ) |>
    
    dplyr::ungroup()
  
  chisq_tab <-
    table(
      soroti_group_composition$Duration_Tercile,
      soroti_group_composition$group_size_cat
    )
  
  list(
    
    group_composition_table =
      group_composition_table,
    
    chisq_table =
      chisq_tab,
    
    chisq_test =
      suppressWarnings(
        chisq.test(chisq_tab)
      )
  )
}


# ------------------------------------------------------------
# Soroti check 2:
# Sample size in each Exposure × Duration tercile cell.
# ------------------------------------------------------------

check_exposure_tenure_cell_sizes <- function(
    TistDat_S
) {
  
  TistDat_S |>
    
    dplyr::mutate(
      
      Exposure_Tercile =
        dplyr::ntile(
          Exposure_sc,
          3
        ),
      
      Exposure_Tercile =
        factor(
          Exposure_Tercile,
          labels = c(
            "Low",
            "Medium",
            "High"
          )
        )
    ) |>
    
    dplyr::count(
      Exposure_Tercile,
      Duration_Tercile
    ) |>
    
    tidyr::pivot_wider(
      names_from = Duration_Tercile,
      values_from = n,
      values_fill = 0
    )
}


# ------------------------------------------------------------
# Translate tenure coefficient into percentage change.
#
# Also returns model-predicted and raw tree counts by tenure
# tercile for checking the "far more trees" claim.
# ------------------------------------------------------------

summarise_tenure_effect_on_trees <- function(
    base_S_trees_nb,
    TistDat_S
) {
  
  tenure_coef_S <-
    summary(base_S_trees_nb)$coefficients$cond[
      "Years_since_reg_sc",
      "Estimate"
    ]
  
  pct_change <-
    100 * (
      exp(tenure_coef_S) - 1
    )
  
  tenure_tercile_means_S <-
    TistDat_S |>
    
    dplyr::group_by(
      Duration_Tercile
    ) |>
    
    dplyr::summarise(
      
      mean_Years_since_reg_sc =
        mean(
          Years_since_reg_sc,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  predicted_trees_by_tenure_S <-
    ggeffects::ggpredict(
      
      base_S_trees_nb,
      
      terms = paste0(
        "Years_since_reg_sc [",
        paste(
          round(
            tenure_tercile_means_S$mean_Years_since_reg_sc,
            3
          ),
          collapse = ","
        ),
        "]"
      )
    )
  
  raw_trees_by_tenure_S <-
    TistDat_S |>
    
    dplyr::group_by(
      Duration_Tercile
    ) |>
    
    dplyr::summarise(
      
      n =
        dplyr::n(),
      
      mean_Trees =
        round(
          mean(
            Trees,
            na.rm = TRUE
          ),
          1
        ),
      
      median_Trees =
        round(
          median(
            Trees,
            na.rm = TRUE
          ),
          1
        ),
      
      .groups = "drop"
    )
  
  list(
    
    tenure_coefficient =
      tenure_coef_S,
    
    pct_change_per_sd =
      pct_change,
    
    tenure_tercile_means =
      tenure_tercile_means_S,
    
    predicted_trees_by_tenure =
      as.data.frame(
        predicted_trees_by_tenure_S
      ),
    
    raw_trees_by_tenure =
      raw_trees_by_tenure_S
  )
}


# ============================================================
# SECTION G: PART 2 -- SIMILARITY MODELS
# H2b / H3b / H4b
# ============================================================


# ------------------------------------------------------------
# Decide which random-effects structure to use for a given
# site/outcome/threshold combination.
#
# Specific evaluation from tested script:
#
# - Density, both sites:
#   always simplified
#
# - Trees, Bushenyi:
#   never simplified
#
# - Trees, Soroti:
#   simplified only at 500m/1000m
# ------------------------------------------------------------

use_simplified_re <- function(
    site,
    outcome,
    threshold_label = NA
) {
  
  if (
    outcome == "Density"
  ) {
    return(TRUE)
  }
  
  if (
    outcome == "Trees" &&
    site == "Bushenyi"
  ) {
    return(FALSE)
  }
  
  if (
    outcome == "Trees" &&
    site == "Soroti"
  ) {
    
    return(
      threshold_label %in%
        c(
          "500m",
          "1000m"
        )
    )
  }
  
  FALSE
}


# ------------------------------------------------------------
# Base similarity model
# ------------------------------------------------------------

fit_base_models_b <- function(
    dat,
    dv_col,
    include_area = FALSE,
    use_simplified = FALSE
) {
  
  re_term <- if (use_simplified) {
    
    "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
    
  } else {
    
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  }
  
  frm <- as.formula(
    paste0(
      dv_col,
      " ~ Exposure_sc + ",
      "Years_since_reg_sc + ",
      "Exposure_sc:Years_since_reg_sc + ",
      "Dist_To_Forest_sc + ",
      if (include_area) {
        "Area_Ha_sc + "
      } else {
        ""
      },
      re_term
    )
  )
  
  lme4::lmer(
    frm,
    data = dat,
    REML = TRUE,
    na.action = na.exclude
  )
}


# ------------------------------------------------------------
# H4b similarity model
# ------------------------------------------------------------

fit_h4b_models <- function(
    dat,
    dv_col,
    nearfar_var,
    include_area = FALSE,
    use_simplified = FALSE
) {
  
  re_term <- if (use_simplified) {
    
    "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
    
  } else {
    
    "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
  }
  
  frm <- as.formula(
    paste0(
      dv_col,
      " ~ Exposure_sc + ",
      "Years_since_reg_sc + ",
      "Exposure_sc:Years_since_reg_sc + ",
      "Dist_To_Forest_sc + ",
      
      if (include_area) {
        "Area_Ha_sc + "
      } else {
        ""
      },
      
      nearfar_var,
      " + Exposure_sc:",
      nearfar_var,
      " + ",
      re_term
    )
  )
  
  lme4::lmer(
    frm,
    data = dat,
    REML = TRUE,
    na.action = na.exclude
  )
}


# ------------------------------------------------------------
# Fit base similarity models across all 10 thresholds.
# One named-list target per site/outcome.
# ------------------------------------------------------------

fit_base_models_b_all_thresholds <- function(
    dat,
    dv_cols,
    site,
    outcome,
    include_area = FALSE,
    threshold_labels = get_threshold_labels()
) {
  
  setNames(
    
    purrr::map2(
      dv_cols,
      threshold_labels,
      
      ~ fit_base_models_b(
        dat,
        dv_col = .x,
        include_area = include_area,
        use_simplified =
          use_simplified_re(
            site,
            outcome,
            .y
          )
      )
    ),
    
    threshold_labels
  )
}


# ------------------------------------------------------------
# Fit H4b similarity models across all 10 thresholds.
# ------------------------------------------------------------

fit_h4b_models_all_thresholds <- function(
    dat,
    dv_cols,
    scaled_nf_cols,
    site,
    outcome,
    include_area = FALSE,
    threshold_labels = get_threshold_labels()
) {
  
  setNames(
    
    purrr::pmap(
      
      list(
        dv_cols,
        scaled_nf_cols,
        threshold_labels
      ),
      
      ~ fit_h4b_models(
        dat,
        dv_col = ..1,
        nearfar_var = ..2,
        include_area = include_area,
        use_simplified =
          use_simplified_re(
            site,
            outcome,
            ..3
          )
      )
    ),
    
    threshold_labels
  )
}


# ------------------------------------------------------------
# Convergence diagnostics for similarity models
# ------------------------------------------------------------

check_convergence_b <- function(
    mod_list,
    label
) {
  
  purrr::map_dfr(
    names(mod_list),
    
    function(th) {
      
      warn <-
        mod_list[[th]]@optinfo$conv$lme4$messages
      
      tibble::tibble(
        
        label = label,
        
        threshold = th,
        
        ok = is.null(warn),
        
        warning =
          ifelse(
            is.null(warn),
            NA_character_,
            paste(
              warn,
              collapse = "; "
            )
          )
      )
    }
  )
}


# ------------------------------------------------------------
# AIC comparison of base-b and H4b models
# ------------------------------------------------------------

compare_aic_b <- function(
    base_list,
    h4b_list,
    threshold_labels,
    site_label
) {
  
  tibble::tibble(
    
    site = site_label,
    
    threshold = threshold_labels,
    
    AIC_base_b =
      purrr::map_dbl(
        base_list,
        AIC
      ),
    
    AIC_h4b =
      purrr::map_dbl(
        h4b_list,
        AIC
      )
    
  ) |>
    
    dplyr::mutate(
      delta_AIC =
        round(
          AIC_base_b - AIC_h4b,
          2
        )
    )
}


# ------------------------------------------------------------
# Extract H2b/H3b coefficients
# ------------------------------------------------------------

extract_h2b_h3b <- function(
    mod_list,
    site_label
) {
  
  purrr::map_dfr(
    names(mod_list),
    
    function(th) {
      
      m <- mod_list[[th]]
      
      cf <-
        as.data.frame(
          coef(summary(m))
        )
      
      cf$term <-
        rownames(cf)
      
      names(cf)[
        names(cf) == "Estimate"
      ] <- "estimate"
      
      names(cf)[
        names(cf) == "Std. Error"
      ] <- "std.error"
      
      names(cf)[
        names(cf) == "t value"
      ] <- "statistic"
      
      cf |>
        
        tibble::as_tibble() |>
        
        dplyr::filter(
          term %in% c(
            "Exposure_sc",
            "Exposure_sc:Years_since_reg_sc"
          )
        ) |>
        
        dplyr::select(
          term,
          estimate,
          std.error,
          statistic
        ) |>
        
        dplyr::mutate(
          threshold = th,
          site = site_label
        )
    }
  )
}


# ------------------------------------------------------------
# Add log1p-transformed tree dissimilarity columns.
# ------------------------------------------------------------

add_log_dissimilarity_columns <- function(
    dat,
    k_range = 1:10
) {
  
  for (k in k_range) {
    
    raw_col <-
      paste0(
        "AbsDissimilarity_Trees_",
        k
      )
    
    log_col <-
      paste0(
        "log_AbsDissimilarity_Trees_",
        k
      )
    
    dat[[log_col]] <-
      log1p(
        dat[[raw_col]]
      )
  }
  
  dat
}


# ------------------------------------------------------------
# Extract Near-Far main effect and interaction from H4b models
# ------------------------------------------------------------

extract_h4b <- function(
    mod_list,
    scaled_nf_cols,
    site_label
) {
  
  purrr::map2_dfr(
    
    names(mod_list),
    scaled_nf_cols,
    
    function(
    th,
    nf_var
    ) {
      
      m <- mod_list[[th]]
      
      cf <-
        as.data.frame(
          coef(summary(m))
        )
      
      cf$term <-
        rownames(cf)
      
      names(cf)[
        names(cf) == "Estimate"
      ] <- "estimate"
      
      names(cf)[
        names(cf) == "Std. Error"
      ] <- "std.error"
      
      names(cf)[
        names(cf) == "t value"
      ] <- "statistic"
      
      interaction_term <-
        paste0(
          "Exposure_sc:",
          nf_var
        )
      
      cf |>
        
        tibble::as_tibble() |>
        
        dplyr::filter(
          term %in% c(
            nf_var,
            interaction_term
          )
        ) |>
        
        dplyr::mutate(
          
          term =
            dplyr::case_when(
              
              term == nf_var ~
                "NearFar_resid_sc",
              
              term == interaction_term ~
                "Exposure_sc:NearFar_resid_sc",
              
              TRUE ~ term
            )
        ) |>
        
        dplyr::select(
          term,
          estimate,
          std.error,
          statistic
        ) |>
        
        dplyr::mutate(
          threshold = th,
          site = site_label
        )
    }
  )
}


# ------------------------------------------------------------
# Fixed-effects covariance between Exposure and Near-Far.
#
# This is NOT raw-data correlation.
# It is the correlation in the fixed-effects covariance matrix.
# ------------------------------------------------------------

extract_h4b_collinearity <- function(
    mod_list,
    scaled_nf_cols,
    site_label
) {
  
  purrr::map2_dfr(
    
    names(mod_list),
    scaled_nf_cols,
    
    function(
    th,
    nf_var
    ) {
      
      m <- mod_list[[th]]
      
      corr_mat <-
        cov2cor(
          as.matrix(
            vcov(m)
          )
        )
      
      tibble::tibble(
        
        threshold = th,
        
        site = site_label,
        
        exposure_nearfar_corr =
          round(
            corr_mat[
              "Exposure_sc",
              nf_var
            ],
            3
          )
      )
    }
  )
}


# ------------------------------------------------------------
# Prepare coefficient table
# ------------------------------------------------------------

prep_coef_table <- function(
    tbl,
    outcome_label,
    threshold_labels = get_threshold_labels()
) {
  
  tbl |>
    
    dplyr::mutate(
      
      Outcome = outcome_label,
      
      threshold =
        factor(
          threshold,
          levels = threshold_labels
        ),
      
      ci_low =
        estimate -
        1.96 * std.error,
      
      ci_high =
        estimate +
        1.96 * std.error
    )
}


# ------------------------------------------------------------
# Plot coefficient estimates by distance threshold
# ------------------------------------------------------------

plot_coef_by_threshold <- function(
    coef_data,
    facet_term = FALSE
) {
  
  p <-
    ggplot2::ggplot(
      coef_data,
      ggplot2::aes(
        x = threshold,
        y = estimate,
        colour = site
      )
    ) +
    
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    ggplot2::geom_pointrange(
      
      ggplot2::aes(
        ymin = ci_low,
        ymax = ci_high
      ),
      
      position =
        ggplot2::position_dodge(
          width = 0.4
        )
    ) +
    
    ggplot2::labs(
      
      x = "Distance threshold",
      
      y = "Coefficient estimate (95% CI)",
      
      colour = "Site"
    ) +
    
    ggplot2::theme_classic(
      base_size = 15
    ) +
    
    ggplot2::theme(
      
      legend.position = "bottom",
      
      legend.text =
        ggplot2::element_text(
          size = 13
        ),
      
      legend.title =
        ggplot2::element_text(
          size = 13
        ),
      
      axis.text.x =
        ggplot2::element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 12
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size = 12
        ),
      
      axis.title =
        ggplot2::element_text(
          size = 14
        ),
      
      strip.text =
        ggplot2::element_text(
          face = "bold",
          size = 13
        )
    )
  
  if (facet_term) {
    
    p +
      ggplot2::facet_grid(
        Outcome ~ term,
        scales = "free_y"
      )
    
  } else {
    
    p +
      ggplot2::facet_wrap(
        ~Outcome,
        scales = "free_y"
      )
  }
}


# ------------------------------------------------------------
# Combine H2b/H3b dissimilarity figures
# ------------------------------------------------------------

combine_h2b_h3b_dissimilarity_figures <- function(
    fig_h2b_exposure,
    fig_h3b_interaction
) {
  
  (
    fig_h2b_exposure /
      fig_h3b_interaction
  ) +
    
    patchwork::plot_layout(
      guides = "collect"
    ) &
    
    ggplot2::theme(
      legend.position = "bottom"
    )
}


# ------------------------------------------------------------
# Norm-convergence asymmetry test.
#
# Fits the same specification separately for farmers above vs.
# below their neighbourhood mean using signed
# Dissimilarity_Trees_k.
#
# This distinguishes genuine two-way convergence from
# one-way "levelling".
# ------------------------------------------------------------

test_convergence_asymmetry <- function(
    dat,
    k,
    site_label,
    threshold_labels = get_threshold_labels(),
    use_simplified = FALSE
) {
  
  signed_col <-
    paste0(
      "Dissimilarity_Trees_",
      k
    )
  
  dat_split <-
    dat |>
    
    dplyr::mutate(
      Above_Neighbour_Mean =
        .data[[signed_col]] > 0
    )
  
  re_term <-
    if (use_simplified) {
      
      "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
      
    } else {
      
      "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
    }
  
  frm <-
    as.formula(
      paste0(
        
        signed_col,
        " ~ Exposure_sc + ",
        "Years_since_reg_sc + ",
        "Exposure_sc:Years_since_reg_sc + ",
        "Dist_To_Forest_sc + ",
        "Area_Ha_sc + ",
        re_term
      )
    )
  
  fit_group <-
    function(
    group_label,
    group_filter
    ) {
      
      sub <-
        dat_split |>
        dplyr::filter(
          Above_Neighbour_Mean ==
            group_filter
        )
      
      if (
        nrow(sub) < 30
      ) {
        return(NULL)
      }
      
      m <-
        tryCatch(
          
          lme4::lmer(
            frm,
            data = sub,
            REML = TRUE,
            na.action = na.exclude
          ),
          
          error = function(e) {
            NULL
          }
        )
      
      if (
        is.null(m)
      ) {
        return(NULL)
      }
      
      cf <-
        as.data.frame(
          coef(summary(m))
        )
      
      cf$term <-
        rownames(cf)
      
      tibble::tibble(
        
        site = site_label,
        
        threshold =
          threshold_labels[k],
        
        group = group_label,
        
        n = nrow(sub),
        
        exposure_est =
          cf$Estimate[
            cf$term ==
              "Exposure_sc"
          ],
        
        exposure_se =
          cf$`Std. Error`[
            cf$term ==
              "Exposure_sc"
          ],
        
        exposure_t =
          cf$`t value`[
            cf$term ==
              "Exposure_sc"
          ]
      )
    }
  
  dplyr::bind_rows(
    
    fit_group(
      "Above neighbour mean",
      TRUE
    ),
    
    fit_group(
      "Below neighbour mean",
      FALSE
    )
  )
}


# ============================================================
# SECTION H: VARIOGRAMS
# STAGE A RAW
# STAGE B RESIDUAL
# ============================================================


# ------------------------------------------------------------
# Generic semivariogram function.
#
# Fits a spherical model.
# Retries with a different starting range if necessary.
# Falls back to a pure nugget if the spherical model cannot
# be reliably estimated.
# ------------------------------------------------------------

fit_variogram_generic <- function(
    data_with_coords,
    value_col,
    cutoff = 8000,
    width = 500
) {
  
  data_clean <-
    data_with_coords |>
    
    dplyr::filter(
      !is.na(
        .data[[value_col]]
      )
    )
  
  data_sf <-
    sf::st_as_sf(
      data_clean,
      coords = c(
        "longitude",
        "latitude"
      ),
      crs = 4326
    )
  
  data_sf <-
    sf::st_transform(
      data_sf,
      32636
    )
  
  data_sp <-
    as(
      data_sf,
      "Spatial"
    )
  
  vgm_empirical <-
    gstat::variogram(
      
      as.formula(
        paste0(
          value_col,
          " ~ 1"
        )
      ),
      
      data_sp,
      
      cutoff = cutoff,
      
      width = width
    )
  
  fit_attempt <-
    function(start_range) {
      
      tryCatch(
        
        gstat::fit.variogram(
          
          vgm_empirical,
          
          gstat::vgm(
            
            psill =
              var(
                data_sp[[value_col]],
                na.rm = TRUE
              ),
            
            model = "Sph",
            
            range = start_range,
            
            nugget = 0
          )
        ),
        
        error = function(e) {
          NULL
        }
      )
    }
  
  vgm_fitted <-
    fit_attempt(
      2000
    )
  
  range_ok <-
    function(v) {
      
      !is.null(v) &&
        
        any(
          v$model == "Sph"
        ) &&
        
        is.finite(
          v$range[
            v$model == "Sph"
          ]
        ) &&
        
        v$range[
          v$model == "Sph"
        ] > 0
    }
  
  if (
    !range_ok(vgm_fitted)
  ) {
    
    vgm_fitted <-
      fit_attempt(
        max(
          vgm_empirical$dist
        ) / 2
      )
  }
  
  if (
    !range_ok(vgm_fitted)
  ) {
    
    vgm_fitted <-
      gstat::vgm(
        
        psill =
          var(
            data_sp[[value_col]],
            na.rm = TRUE
          ),
        
        model = "Nug"
      )
  }
  
  fit_curve <-
    gstat::variogramLine(
      vgm_fitted,
      maxdist =
        max(
          vgm_empirical$dist
        )
    )
  
  nugget_val <-
    if (
      any(
        vgm_fitted$model == "Nug"
      )
    ) {
      
      vgm_fitted$psill[
        vgm_fitted$model == "Nug"
      ]
      
    } else {
      
      0
    }
  
  partial_sill_val <-
    if (
      any(
        vgm_fitted$model == "Sph"
      )
    ) {
      
      vgm_fitted$psill[
        vgm_fitted$model == "Sph"
      ]
      
    } else {
      
      0
    }
  
  range_val <-
    if (
      any(
        vgm_fitted$model == "Sph"
      )
    ) {
      
      vgm_fitted$range[
        vgm_fitted$model == "Sph"
      ]
      
    } else {
      
      0
    }
  
  total_sill_val <-
    nugget_val +
    partial_sill_val
  
  spatial_dependence_val <-
    if (
      total_sill_val > 0
    ) {
      
      partial_sill_val /
        total_sill_val
      
    } else {
      
      0
    }
  
  list(
    
    empirical =
      vgm_empirical,
    
    fitted =
      vgm_fitted,
    
    fit_curve =
      fit_curve,
    
    nugget =
      nugget_val,
    
    partial_sill =
      partial_sill_val,
    
    total_sill =
      total_sill_val,
    
    range_m =
      range_val,
    
    spatial_dependence =
      spatial_dependence_val
  )
}


# ------------------------------------------------------------
# Draw one semivariogram panel on an existing base-R device.
# ------------------------------------------------------------

draw_variogram_panel <- function(
    vgm_result,
    plot_title,
    label = NULL,
    show_stats = FALSE
) {
  
  dist_full <-
    c(
      0,
      vgm_result$empirical$dist
    )
  
  gamma_full <-
    c(
      vgm_result$nugget,
      vgm_result$empirical$gamma
    )
  
  y_max <-
    max(
      gamma_full,
      vgm_result$total_sill * 1.1
    )
  
  plot(
    
    dist_full,
    gamma_full,
    
    type = "b",
    
    pch = 19,
    
    col = "black",
    
    xlim =
      c(
        0,
        max(dist_full)
      ),
    
    ylim =
      c(
        0,
        y_max
      ),
    
    xlab =
      "Distance (m)",
    
    ylab =
      "Semivariance",
    
    main =
      plot_title
  )
  
  lines(
    
    vgm_result$fit_curve$dist,
    
    vgm_result$fit_curve$gamma,
    
    col = "gray40",
    
    lwd = 2
  )
  
  abline(
    
    v = vgm_result$range_m,
    
    lty = 2,
    
    lwd = 2,
    
    col = "red"
  )
  
  text(
    
    x = vgm_result$range_m,
    
    y = y_max * 0.85,
    
    labels =
      paste0(
        "Range = ",
        round(
          vgm_result$range_m,
          0
        ),
        " m"
      ),
    
    col = "red",
    
    cex = 1.1,
    
    pos = 4
  )
  
  if (
    show_stats
  ) {
    
    spatial_pct <-
      round(
        vgm_result$spatial_dependence * 100,
        0
      )
    
    unexplained_pct <-
      100 -
      spatial_pct
    
    stat_x <-
      par("usr")[2] -
      0.03 *
      diff(
        par("usr")[1:2]
      )
    
    stat_y <-
      y_max * 0.15
    
    text(
      
      x = stat_x,
      
      y = stat_y,
      
      labels =
        paste0(
          "Spatial dependence = ",
          spatial_pct,
          "%\n(",
          unexplained_pct,
          "% unexplained)"
        ),
      
      col = "black",
      
      cex = 1.0,
      
      adj = c(
        1,
        0
      )
    )
  }
  
  if (
    !is.null(label)
  ) {
    
    text(
      
      x =
        par("usr")[1],
      
      y =
        par("usr")[4],
      
      labels =
        label,
      
      font = 2,
      
      cex = 1.6,
      
      adj =
        c(
          -0.3,
          1.3
        ),
      
      xpd = NA
    )
  }
}


# ------------------------------------------------------------
# Standard single variogram plot
# ------------------------------------------------------------

plot_variogram_standard <- function(
    vgm_result,
    plot_title,
    save_path,
    show_stats = FALSE
) {
  
  dir.create(
    dirname(save_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  png(
    save_path,
    width = 2400,
    height = 1800,
    res = 300
  )
  
  par(
    mar = c(
      5,
      5,
      3,
      2
    ),
    cex.axis = 1.3,
    cex.lab = 1.4,
    cex.main = 1.5
  )
  
  draw_variogram_panel(
    vgm_result,
    plot_title,
    show_stats = show_stats
  )
  
  dev.off()
  
  save_path
}


# ------------------------------------------------------------
# Multi-panel variogram figure
# ------------------------------------------------------------

make_variogram_panel_figure <- function(
    results_list,
    titles,
    save_path,
    nrow = 2,
    ncol = 2,
    labels = LETTERS[
      seq_along(results_list)
    ],
    show_stats = FALSE
) {
  
  dir.create(
    dirname(save_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  png(
    
    save_path,
    
    width =
      1600 * ncol,
    
    height =
      1200 * nrow,
    
    res = 300
  )
  
  par(
    
    mfrow =
      c(
        nrow,
        ncol
      ),
    
    mar =
      c(
        5,
        5,
        3,
        2
      ),
    
    cex.axis = 1.1,
    
    cex.lab = 1.2,
    
    cex.main = 1.3
  )
  
  for (
    i in seq_along(results_list)
  ) {
    
    draw_variogram_panel(
      
      results_list[[i]],
      
      titles[i],
      
      label =
        labels[i],
      
      show_stats =
        show_stats
    )
  }
  
  dev.off()
  
  save_path
}


# ------------------------------------------------------------
# Stage A:
# Raw exploratory variogram for one site/outcome.
#
# By default outcome is log1p transformed.
# ------------------------------------------------------------

run_raw_variogram <- function(
    full_data,
    site_name,
    outcome_col,
    outcome_label,
    output_dir,
    log_transform = TRUE,
    show_stats = FALSE
) {
  
  site_data <-
    full_data |>
    dplyr::filter(
      Proj_Area == site_name
    )
  
  value_col <-
    outcome_col
  
  if (
    log_transform
  ) {
    
    value_col <-
      paste0(
        "log_",
        outcome_col
      )
    
    site_data[[value_col]] <-
      log1p(
        site_data[[outcome_col]]
      )
  }
  
  vgm_result <-
    fit_variogram_generic(
      site_data,
      value_col
    )
  
  plot_title <-
    paste0(
      site_name,
      " - ",
      outcome_label
    )
  
  save_path <-
    file.path(
      
      output_dir,
      
      paste0(
        site_name,
        "_variogram_",
        outcome_col,
        ".png"
      )
    )
  
  plot_variogram_standard(
    
    vgm_result,
    
    plot_title,
    
    save_path,
    
    show_stats =
      show_stats
  )
  
  vgm_result
}


# ------------------------------------------------------------
# Raw variogram summary
# ------------------------------------------------------------

build_raw_variogram_summary <- function(
    density_B,
    density_S,
    trees_B,
    trees_S
) {
  
  tibble::tibble(
    
    site = c(
      "Bushenyi",
      "Soroti",
      "Bushenyi",
      "Soroti"
    ),
    
    outcome = c(
      "Density (log)",
      "Density (log)",
      "Trees (log)",
      "Trees (log)"
    ),
    
    nugget = c(
      density_B$nugget,
      density_S$nugget,
      trees_B$nugget,
      trees_S$nugget
    ),
    
    total_sill = c(
      density_B$total_sill,
      density_S$total_sill,
      trees_B$total_sill,
      trees_S$total_sill
    ),
    
    range_m = c(
      density_B$range_m,
      density_S$range_m,
      trees_B$range_m,
      trees_S$range_m
    ),
    
    spatial_dependence_pct =
      round(
        100 *
          c(
            density_B$spatial_dependence,
            density_S$spatial_dependence,
            trees_B$spatial_dependence,
            trees_S$spatial_dependence
          ),
        0
      )
  )
}


# ------------------------------------------------------------
# Stage B:
# Exposure-only model.
#
# Fixed effect:
#   Exposure_sc
#
# Same RE hierarchy as elsewhere.
# ------------------------------------------------------------

fit_exposure_only_model <- function(
    dat,
    dv_col,
    use_simplified = FALSE
) {
  
  re_term <-
    if (
      use_simplified
    ) {
      
      "(1 | Village_ID) + (1 | Cluster_ID/Group_ID)"
      
    } else {
      
      "(1 | Admin_Districts/Subcounty/Village_ID) + (1 | Cluster_ID/Group_ID)"
    }
  
  frm <-
    as.formula(
      paste0(
        dv_col,
        " ~ Exposure_sc + ",
        re_term
      )
    )
  
  lme4::lmer(
    frm,
    data = dat,
    REML = TRUE,
    na.action = na.exclude
  )
}


# ------------------------------------------------------------
# Exposure residual variogram
# ------------------------------------------------------------

run_exposure_residual_variogram <- function(
    dat,
    outcome_col,
    site_name,
    outcome_label,
    use_simplified = FALSE
) {
  
  dat_complete <-
    dat |>
    
    dplyr::filter(
      
      !is.na(longitude),
      
      !is.na(latitude),
      
      !is.na(
        .data[[outcome_col]]
      )
    )
  
  log_col <-
    paste0(
      "log_",
      outcome_col
    )
  
  dat_complete[[log_col]] <-
    log1p(
      dat_complete[[outcome_col]]
    )
  
  m_exposure <-
    fit_exposure_only_model(
      
      dat_complete,
      
      log_col,
      
      use_simplified =
        use_simplified
    )
  
  dat_res <-
    dat_complete |>
    
    dplyr::mutate(
      resid_exposure =
        resid(m_exposure)
    )
  
  fit_variogram_generic(
    
    dat_res |>
      dplyr::rename(
        resid_val =
          resid_exposure
      ),
    
    "resid_val"
  )
}


# ------------------------------------------------------------
# Exposure residual variogram summary
# ------------------------------------------------------------

build_exposure_residual_summary <- function(
    density_B,
    density_S,
    trees_B,
    trees_S
) {
  
  tibble::tibble(
    
    site = c(
      "Bushenyi",
      "Soroti",
      "Bushenyi",
      "Soroti"
    ),
    
    outcome = c(
      "Density (log)",
      "Density (log)",
      "Trees (log)",
      "Trees (log)"
    ),
    
    nugget = c(
      density_B$nugget,
      density_S$nugget,
      trees_B$nugget,
      trees_S$nugget
    ),
    
    total_sill = c(
      density_B$total_sill,
      density_S$total_sill,
      trees_B$total_sill,
      trees_S$total_sill
    ),
    
    range_m = c(
      density_B$range_m,
      density_S$range_m,
      trees_B$range_m,
      trees_S$range_m
    ),
    
    spatial_dependence_pct =
      round(
        100 *
          c(
            density_B$spatial_dependence,
            density_S$spatial_dependence,
            trees_B$spatial_dependence,
            trees_S$spatial_dependence
          ),
        0
      )
  )
}


# ------------------------------------------------------------
# Variance explained from reduction in semivariogram sill.
#
# Formula:
#
#   (sill_raw - sill_residual) / sill_raw
#
# Tables are joined by row position, matching the tested
# script's bind_cols() approach.
# ------------------------------------------------------------

build_variance_explained_table <- function(
    raw_summary,
    residual_summary
) {
  
  raw_summary |>
    
    dplyr::select(
      site,
      outcome,
      sill_raw = total_sill
    ) |>
    
    dplyr::bind_cols(
      
      residual_summary |>
        dplyr::select(
          sill_residual =
            total_sill
        )
    ) |>
    
    dplyr::mutate(
      
      pct_variance_explained_by_exposure =
        round(
          100 *
            (
              sill_raw -
                sill_residual
            ) /
            sill_raw,
          1
        )
    )
}


# ------------------------------------------------------------
# Word-ready residual comparison table.
#
# Uses flextable/officer-compatible output.
# ------------------------------------------------------------

build_residual_comparison_docx <- function(
    raw_summary,
    residual_summary,
    variance_explained,
    range_resolved,
    output_path
) {
  
  residual_comparison_table <-
    raw_summary |>
    
    dplyr::select(
      
      site,
      
      outcome,
      
      range_raw =
        range_m,
      
      spatial_dependence_raw =
        spatial_dependence_pct
    ) |>
    
    dplyr::bind_cols(
      
      residual_summary |>
        dplyr::select(
          
          range_residual =
            range_m,
          
          spatial_dependence_residual =
            spatial_dependence_pct
        )
    ) |>
    
    dplyr::bind_cols(
      
      variance_explained |>
        dplyr::select(
          pct_variance_explained_by_exposure
        )
    ) |>
    
    dplyr::mutate(
      
      range_residual_display =
        ifelse(
          
          range_resolved,
          
          paste0(
            round(
              range_residual,
              0
            ),
            " m"
          ),
          
          "Not resolved†"
        ),
      
      range_raw_display =
        paste0(
          round(
            range_raw,
            0
          ),
          " m"
        ),
      
      outcome =
        stringr::str_remove(
          outcome,
          " \\(log\\)"
        )
    ) |>
    
    dplyr::select(
      
      Site = site,
      
      Outcome = outcome,
      
      `Range, raw` =
        range_raw_display,
      
      `Range, residual` =
        range_residual_display,
      
      `Spatial dependence, raw (%)` =
        spatial_dependence_raw,
      
      `Spatial dependence, residual (%)` =
        spatial_dependence_residual,
      
      `Variance explained by exposure (%)` =
        pct_variance_explained_by_exposure
    )
  
  ft <-
    flextable::flextable(
      residual_comparison_table
    ) |>
    
    flextable::set_caption(
      "Table SX. Comparison of raw and exposure-only residual semivariogram parameters, by site and outcome."
    ) |>
    
    flextable::add_footer_lines(
      
      "† Fitted spherical range exceeded the maximum tested distance (8,000 m); the empirical semivariogram had not plateaued within the observed range. Spatial dependence values for these rows should be interpreted with caution; variance explained by exposure (final column) does not share this limitation, as it is computed from total sills rather than the plateau location."
    ) |>
    
    flextable::theme_vanilla() |>
    
    flextable::autofit() |>
    
    flextable::fontsize(
      size = 9,
      part = "all"
    ) |>
    
    flextable::bold(
      part = "header"
    )
  
  dir.create(
    dirname(output_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  flextable::save_as_docx(
    ft,
    path = output_path
  )
  
  output_path
}


# ============================================================
# SECTION I: PART 5 -- SPATIAL DIAGNOSTICS
# MORAN'S I / DHARMa
# ============================================================


# ------------------------------------------------------------
# Extract residuals from density and tree models
# ------------------------------------------------------------

extract_model_residuals <- function(
    dat,
    density_model,
    trees_model
) {
  
  dat$res_density <-
    resid(
      density_model
    )
  
  dat$res_trees <-
    residuals(
      trees_model,
      type = "pearson"
    )
  
  dat
}


# ------------------------------------------------------------
# Aggregate model residuals to village level
# ------------------------------------------------------------

aggregate_residuals_to_village <- function(
    dat_with_resid
) {
  
  dat_with_resid |>
    
    dplyr::group_by(
      Village_ID
    ) |>
    
    dplyr::summarise(
      
      res_density =
        mean(
          res_density,
          na.rm = TRUE
        ),
      
      res_trees =
        mean(
          res_trees,
          na.rm = TRUE
        ),
      
      lon =
        mean(
          longitude,
          na.rm = TRUE
        ),
      
      lat =
        mean(
          latitude,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
}


# ------------------------------------------------------------
# Build k-nearest-neighbour spatial weights.
#
# k = 4 by default because villages are irregularly spaced
# rather than sharing contiguous borders.
# ------------------------------------------------------------

build_spatial_weights <- function(
    village_res,
    utm_crs = 32636,
    k = 4
) {
  
  sf_obj <-
    sf::st_as_sf(
      
      village_res,
      
      coords =
        c(
          "lon",
          "lat"
        ),
      
      crs = 4326
    )
  
  sf_proj <-
    sf::st_transform(
      sf_obj,
      utm_crs
    )
  
  coords <-
    sf::st_coordinates(
      sf_proj
    )
  
  nb <-
    spdep::knn2nb(
      spdep::knearneigh(
        coords,
        k = k
      )
    )
  
  lw <-
    spdep::nb2listw(
      nb,
      style = "W"
    )
  
  list(
    
    sf =
      sf_proj,
    
    coords =
      coords,
    
    nb =
      nb,
    
    lw =
      lw
  )
}


# ------------------------------------------------------------
# Moran's I
# ------------------------------------------------------------

run_morans_i <- function(
    village_res,
    lw,
    value_col
) {
  
  spdep::moran.test(
    village_res[[value_col]],
    lw
  )
}


# ------------------------------------------------------------
# DHARMa simulated residuals
# ------------------------------------------------------------

run_dharma_diagnostics <- function(
    model,
    n = 1000
) {
  
  DHARMa::simulateResiduals(
    model,
    n = n
  )
}


# ------------------------------------------------------------
# DHARMa spatial autocorrelation test
# ------------------------------------------------------------

run_spatial_autocorrelation_test <- function(
    sim,
    x,
    y,
    plot = TRUE
) {
  
  DHARMa::testSpatialAutocorrelation(
    sim,
    x = x,
    y = y,
    plot = plot
  )
}


# ------------------------------------------------------------
# DHARMa dispersion test
# ------------------------------------------------------------

run_dispersion_test <- function(
    sim
) {
  
  DHARMa::testDispersion(
    sim
  )
}


# ------------------------------------------------------------
# DHARMa zero-inflation test
# ------------------------------------------------------------

run_zeroinflation_test <- function(
    sim
) {
  
  DHARMa::testZeroInflation(
    sim
  )
}