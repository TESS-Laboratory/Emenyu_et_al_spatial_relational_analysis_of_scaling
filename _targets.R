# Created by use_targets().
# Follow the comments below to fill in each target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(magrittr)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "ergm.multi", "network",
               "geosphere", "ergm", "sna", "statnet","lme4","glmmTMB","broom.mixed",
               "patchwork", "modelsummary", "terra","sf", "scales","ggridges" )
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Created by use_targets().
#
# Complete pipeline:
#   - TIST data
#   - ESA WorldCover landcover
#   - landcover suitability
#   - eligible land / forest exclusion
#   - project areas
#   - TIST/random points
#   - forest-distance analysis
#   - ERGM network analysis
#   - ERGM diagnostics
#   - overlap simulations
#   - publication-ready tables
#   - mixed-effects outcome models


################################################################################
# SOURCE FUNCTIONS
################################################################################

tar_source()


list(
  
  
  ################################################################################
  # DATA
  ################################################################################
  
  tar_target(
    tist_data,
    load_TIST_data(
      "Data/TISTDat/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv"
    )
  ),
  
  
  ################################################################################
  # STATIC SPATIAL DATA
  ################################################################################
  
  tar_target(
    Uganda,
    load_uganda_boundary(
      "Data/UGshapefiles/uga_admbnda_adm0_ubos_20200824.shp"
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # ESA WorldCover 2020 tiles covering Uganda
  # --------------------------------------------------------------------------
  
  tar_target(
    landcover_files,
    c(
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
  ),
  
  
  # --------------------------------------------------------------------------
  # Persistent land-cover mosaic on disk
  # --------------------------------------------------------------------------
  
  tar_target(
    landcover_raster_file,
    build_landcover_mosaic_file(
      files = landcover_files,
      boundary = Uganda,
      output =
        "Data/landcover/ESA_WorldCover_2020/Uganda_WorldCover_2020_mosaic.tif"
    ),
    format = "file"
  ),
  
  
  tar_target(
    landcover_legend,
    create_landcover_legend()
  ),
  
  
  ################################################################################
  # LAND-COVER PROCESSING
  ################################################################################
  
  tar_target(
    TIST_landcover,
    extract_landcover(
      tist_data,
      landcover_raster_file
    )
  ),
  
  
  tar_target(
    Bushenyi_landcover,
    summarise_landcover(
      TIST_landcover,
      "Bushenyi",
      landcover_legend
    )
  ),
  
  
  tar_target(
    Soroti_landcover,
    summarise_landcover(
      TIST_landcover,
      "Soroti",
      landcover_legend
    )
  ),
  
  
  tar_target(
    landcover_typology,
    create_landcover_typology()
  ),
  
  
  tar_target(
    Bushenyi_suitability,
    apply_tree_suitability(
      Bushenyi_landcover,
      landcover_typology
    )
  ),
  
  
  tar_target(
    Soroti_suitability,
    apply_tree_suitability(
      Soroti_landcover,
      landcover_typology
    )
  ),
  
  
  ################################################################################
  # ELIGIBLE LAND / FOREST EXCLUSION
  ################################################################################
  
  tar_target(
    eligible_mask_file,
    create_landcover_mask_file(
      landcover_raster_file
    ),
    format = "file"
  ),
  
  
  tar_target(
    forest_reserves,
    sf::st_read(
      "Data/Uganda_Forest_Reserves/Uganda_Forest_Reserves.shp",
      quiet = TRUE
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Forest reserves prepared as a persistent vector file
  # in the land-cover raster CRS for rasterisation
  # --------------------------------------------------------------------------
  
  tar_target(
    forest_reserves_vector_file,
    prepare_forest_reserves_for_exclusion_file(
      forest_reserves,
      landcover_raster_file
    ),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Forest reserves prepared as a persistent vector file
  # in UTM Zone 36N for distance calculations
  # --------------------------------------------------------------------------
  
  tar_target(
    forest_reserves_projected_file,
    prepare_forest_reserves_file(
      forest_reserves,
      target_crs = 32636
    ),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Rasterise forest reserves
  # --------------------------------------------------------------------------
  
  tar_target(
    forest_mask_file,
    create_forest_exclusion_mask_file(
      forest_reserves_vector_file,
      landcover_raster_file
    ),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Remove forest reserves from eligible land
  # --------------------------------------------------------------------------
  
  tar_target(
    eligible_area_file,
    exclude_forest_reserves_file(
      eligible_mask_file,
      forest_mask_file
    ),
    format = "file"
  ),
  
  
  ################################################################################
  # PROJECT AREAS
  ################################################################################
  
  # --------------------------------------------------------------------------
  # Uganda ADM2 boundaries
  # --------------------------------------------------------------------------
  
  tar_target(
    districts_sf,
    sf::st_read(
      "Data/UGshapefiles/uga_admbnda_adm2_ubos_20200824.shp",
      quiet = TRUE
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Bushenyi district names
  #
  # These MUST match ADM2_EN in districts_sf.
  # --------------------------------------------------------------------------
  
  tar_target(
    Bushenyi_names,
    c(
      "Buhweju",
      "Bushenyi",
      "Ibanda",
      "Kamwenge",
      "Kasese",
      "Kitagwenda",
      "Kyenjojo",
      "Mitooma",
      "Ntungamo",
      "Rubirizi",
      "Rwampara",
      "Sheema"
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Soroti district names
  # --------------------------------------------------------------------------
  
  tar_target(
    Soroti_names,
    c(
      "Alebtong",
      "Amuria",
      "Kalaki",
      "Kapelebyong",
      "Serere",
      "Soroti"
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Bushenyi project area
  # --------------------------------------------------------------------------
  
  tar_target(
    Bushenyi_area,
    create_project_area(
      districts_sf,
      Bushenyi_names
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Soroti project area
  # --------------------------------------------------------------------------
  
  tar_target(
    Soroti_area,
    create_project_area(
      districts_sf,
      Soroti_names
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Eligible Bushenyi area
  # --------------------------------------------------------------------------
  
  tar_target(
    eligible_Bushenyi,
    mask_project_area_file(
      eligible_area_file,
      Bushenyi_area,
      output = "Data/derived/eligible_bushenyi.tif"
    ),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Eligible Soroti area
  # --------------------------------------------------------------------------
  
  tar_target(
    eligible_Soroti,
    mask_project_area_file(
      eligible_area_file,
      Soroti_area,
      output = "Data/derived/eligible_soroti.tif"
    ),
    format = "file"
  ),
  
  
  ################################################################################
  # FOREST-DISTANCE ANALYSIS
  ################################################################################
  
  tar_target(
    TIST_Bushenyi_points,
    prepare_TIST_distance_points(
      tist_data,
      "Bushenyi"
    )
  ),
  
  
  tar_target(
    TIST_Soroti_points,
    prepare_TIST_distance_points(
      tist_data,
      "Soroti"
    )
  ),
  
  
  tar_target(
    random_Bushenyi_points,
    prepare_random_ag_points(
      eligible_Bushenyi,
      n_points = 1000
    )
  ),
  
  
  tar_target(
    random_Soroti_points,
    prepare_random_ag_points(
      eligible_Soroti,
      n_points = 1000
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Bushenyi forest distances
  # --------------------------------------------------------------------------
  
  tar_target(
    Bushenyi_TIST_distance,
    nearest_forest_distance(
      TIST_Bushenyi_points,
      forest_reserves_projected_file
    )
  ),
  
  
  tar_target(
    Bushenyi_random_distance,
    nearest_forest_distance(
      random_Bushenyi_points,
      forest_reserves_projected_file
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Soroti forest distances
  # --------------------------------------------------------------------------
  
  tar_target(
    Soroti_TIST_distance,
    nearest_forest_distance(
      TIST_Soroti_points,
      forest_reserves_projected_file
    )
  ),
  
  
  tar_target(
    Soroti_random_distance,
    nearest_forest_distance(
      random_Soroti_points,
      forest_reserves_projected_file
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Distance comparisons
  # --------------------------------------------------------------------------
  
  tar_target(
    Bushenyi_distance_comparison,
    create_distance_comparison(
      Bushenyi_TIST_distance,
      Bushenyi_random_distance
    )
  ),
  
  
  tar_target(
    Soroti_distance_comparison,
    create_distance_comparison(
      Soroti_TIST_distance,
      Soroti_random_distance
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Distance summaries
  # --------------------------------------------------------------------------
  
  tar_target(
    Bushenyi_distance_summary,
    summarise_forest_distance(
      Bushenyi_distance_comparison
    )
  ),
  
  
  tar_target(
    Soroti_distance_summary,
    summarise_forest_distance(
      Soroti_distance_comparison
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Distance tests
  # --------------------------------------------------------------------------
  
  tar_target(
    Bushenyi_distance_test,
    test_forest_distance(
      Bushenyi_distance_comparison
    )
  ),
  
  
  tar_target(
    Soroti_distance_test,
    test_forest_distance(
      Soroti_distance_comparison
    )
  ),
  
  
  ################################################################################
  # LAND-COVER FIGURES
  ################################################################################
  
  tar_target(
    landcover_colors,
    get_landcover_colors()
  ),
  
  
  tar_target(
    Bushenyi_landcover_plot,
    plot_landcover_pie(
      Bushenyi_landcover,
      "Bushenyi",
      landcover_colors
    )
  ),
  
  
  tar_target(
    Soroti_landcover_plot,
    plot_landcover_pie(
      Soroti_landcover,
      "Soroti",
      landcover_colors
    )
  ),
  
  
  tar_target(
    landcover_comparison_figure,
    combine_two_area_plots(
      Bushenyi_landcover_plot,
      Soroti_landcover_plot
    )
  ),
  
  
  tar_target(
    Bushenyi_suitability_plot,
    plot_suitability_bar(
      Bushenyi_suitability,
      "Bushenyi"
    )
  ),
  
  
  tar_target(
    Soroti_suitability_plot,
    plot_suitability_bar(
      Soroti_suitability,
      "Soroti"
    )
  ),
  
  
  ################################################################################
  # DISTANCE-SPECIFIC NETWORK GRID
  ################################################################################
  
  tar_target(
    sites,
    c(
      "Bushenyi",
      "Soroti"
    )
  ),
  
  
  tar_target(
    distance_thresholds,
    c(
      500,
      2000,
      3500
    )
  ),
  
  
  tar_target(
    site_distance_grid,
    expand.grid(
      site = sites,
      distance_threshold = distance_thresholds
    )
  ),
  
  
  ################################################################################
  # NETWORK GENERATION
  ################################################################################
  
  tar_target(
    site_networks,
    build_site_networks(
      tist_data,
      site = site_distance_grid$site,
      distance_threshold =
        site_distance_grid$distance_threshold
    ),
    pattern = map(site_distance_grid)
  ),
  
  
  ################################################################################
  # ERGM MODELS
  ################################################################################
  
  tar_target(
    site_models,
    fit_site_model(
      site_networks
    ),
    pattern = map(site_networks)
  ),
  
  
  ################################################################################
  # MODEL DIAGNOSTICS
  ################################################################################
  
  tar_target(
    model_diagnostics,
    {
      
      cat(
        "\n==============================\n",
        "Site:",
        site_models$site,
        "\nDistance:",
        site_models$distance_threshold,
        "m\n",
        "==============================\n"
      )
      
      print(
        summary(
          site_models$model
        )
      )
      
      site_models
      
    },
    pattern = map(site_models)
  ),
  
  
  ################################################################################
  # ERGM OVERLAP SIMULATIONS
  ################################################################################
  
  tar_target(
    overlap_simulations,
    evaluate_overlap(
      site_models,
      nsim = 1000
    ),
    pattern = map(site_models)
  ),
  
  
  ################################################################################
  # SIMULATION DIAGNOSTICS
  ################################################################################
  
  tar_target(
    simulation_diagnostics,
    {
      
      cat(
        "\n==============================\n",
        "Site:",
        overlap_simulations$site,
        "\nDistance:",
        overlap_simulations$distance_threshold,
        "m\n",
        "==============================\n"
      )
      
      print(
        overlap_simulations$observed
      )
      
      print(
        overlap_simulations$p.values
      )
      
      overlap_simulations
      
    },
    pattern = map(overlap_simulations)
  ),
  
  
  ################################################################################
  # PUBLICATION-READY ERGM RESULTS
  ################################################################################
  
  tar_target(
    ergm_results,
    {
      
      coef <- summary(
        site_models$model
      )$coefficients
      
      tibble(
        Site = site_models$site,
        Distance = site_models$distance_threshold,
        Term = rownames(coef),
        Estimate = coef[, 1],
        Std_Error = coef[, 2],
        Z_value = coef[, 3],
        P_value = coef[, 4]
      )
      
    },
    pattern = map(site_models)
  ),
  
  
  ################################################################################
  # PUBLICATION-READY OVERLAP RESULTS
  ################################################################################
  
  tar_target(
    overlap_results,
    {
      
      tibble(
        
        Site =
          overlap_simulations$site,
        
        Distance =
          overlap_simulations$distance_threshold,
        
        Statistic =
          c(
            "gcor",
            "Overlap"
          ),
        
        Observed =
          c(
            overlap_simulations$observed$gcor,
            overlap_simulations$observed$overlap
          ),
        
        Simulated_mean =
          c(
            mean(
              overlap_simulations$simulated$gcor
            ),
            mean(
              overlap_simulations$simulated$overlap
            )
          ),
        
        Simulated_SD =
          c(
            sd(
              overlap_simulations$simulated$gcor
            ),
            sd(
              overlap_simulations$simulated$overlap
            )
          ),
        
        P_value =
          c(
            overlap_simulations$p.values$gcor,
            overlap_simulations$p.values$overlap
          )
      )
      
    },
    pattern = map(overlap_simulations)
  ),
  
  
  ################################################################################
  # MIXED-EFFECTS OUTCOME MODELS
  ################################################################################
  
  
  tar_target(
    model_data,
    tist_data %>%
      filter(
        Proj_Area %in% c(
          "Bushenyi",
          "Soroti"
        )
      ) %>%
      mutate(
        Proj_Area = factor(Proj_Area),
        
        Exposure_sc =
          as.numeric(scale(Exposure)),
        
        Years_since_reg_sc =
          as.numeric(scale(Years_since_reg)),
        
        Dist_To_Forest_sc =
          as.numeric(scale(Dist_To_Forest_m)),
        
        across(
          starts_with("NearFar_resid_"),
          ~ as.numeric(scale(.)),
          .names = "{.col}_sc"
        )
      )
  ),
  
  
  tar_target(
    Bushenyi_model_data,
    model_data %>%
      filter(
        Proj_Area == "Bushenyi"
      )
  ),
  
  
  tar_target(
    Soroti_model_data,
    model_data %>%
      filter(
        Proj_Area == "Soroti"
      )
  ),
  
  
  ################################################################################
  # BASE MODELS
  ################################################################################
  
  
  tar_target(
    Bushenyi_base_density_model,
    lme4::lmer(
      Density_winsor99 ~
        Exposure_sc +
        Years_since_reg_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Village_ID) +
        (1 | Cluster_ID) +
        (1 | Group_ID),
      data = Bushenyi_model_data,
      REML = TRUE
    )
  ),
  
  
  tar_target(
    Soroti_base_density_model,
    lme4::lmer(
      Density_winsor99 ~
        Exposure_sc +
        Years_since_reg_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Village_ID) +
        (1 | Cluster_ID) +
        (1 | Group_ID),
      data = Soroti_model_data,
      REML = TRUE
    )
  ),
  
  
  tar_target(
    Bushenyi_base_trees_model,
    glmmTMB::glmmTMB(
      Trees ~
        Exposure_sc +
        Years_since_reg_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Admin_Districts / Subcounty / Village_ID) +
        (1 | Cluster_ID / Group_ID),
      data = Bushenyi_model_data,
      family = glmmTMB::nbinom2
    )
  ),
  
  
  tar_target(
    Soroti_base_trees_model,
    glmmTMB::glmmTMB(
      Trees ~
        Exposure_sc +
        Years_since_reg_sc +
        Exposure_sc:Years_since_reg_sc +
        Dist_To_Forest_sc +
        (1 | Admin_Districts / Subcounty / Village_ID) +
        (1 | Cluster_ID / Group_ID),
      data = Soroti_model_data,
      family = glmmTMB::nbinom2
    )
  ),
  
  
  ################################################################################
  # H4 NEAR-FAR MODELS
  ################################################################################
  
  
  tar_target(
    H4_thresholds,
    seq(
      1,
      10
    )
  ),
  
  
  tar_target(
    Bushenyi_H4_density_models,
    {
      
      purrr::map(
        H4_thresholds,
        function(i) {
          
          predictor <- paste0(
            "NearFar_resid_",
            i,
            "_sc"
          )
          
          interaction_term <- paste0(
            "Exposure_sc:",
            predictor
          )
          
          formula_text <- paste(
            "Density_winsor99 ~",
            "Exposure_sc +",
            "Years_since_reg_sc +",
            "Exposure_sc:Years_since_reg_sc +",
            "Dist_To_Forest_sc +",
            predictor,
            "+",
            interaction_term,
            "+ (1 | Village_ID)",
            "+ (1 | Cluster_ID)",
            "+ (1 | Group_ID)"
          )
          
          lme4::lmer(
            as.formula(formula_text),
            data = Bushenyi_model_data,
            REML = TRUE
          )
          
        }
      )
      
    }
  ),
  
  
  tar_target(
    Soroti_H4_density_models,
    {
      
      purrr::map(
        H4_thresholds,
        function(i) {
          
          predictor <- paste0(
            "NearFar_resid_",
            i,
            "_sc"
          )
          
          interaction_term <- paste0(
            "Exposure_sc:",
            predictor
          )
          
          formula_text <- paste(
            "Density_winsor99 ~",
            "Exposure_sc +",
            "Years_since_reg_sc +",
            "Exposure_sc:Years_since_reg_sc +",
            "Dist_To_Forest_sc +",
            predictor,
            "+",
            interaction_term,
            "+ (1 | Admin_Districts / Subcounty / Village_ID)",
            "+ (1 | Cluster_ID / Group_ID)"
          )
          
          lme4::lmer(
            as.formula(formula_text),
            data = Soroti_model_data,
            REML = TRUE
          )
          
        }
      )
      
    }
  ),
  
  
  tar_target(
    Bushenyi_H4_trees_models,
    {
      
      purrr::map(
        H4_thresholds,
        function(i) {
          
          predictor <- paste0(
            "NearFar_resid_",
            i,
            "_sc"
          )
          
          interaction_term <- paste0(
            "Exposure_sc:",
            predictor
          )
          
          formula_text <- paste(
            "Trees ~",
            "Exposure_sc +",
            "Years_since_reg_sc +",
            "Exposure_sc:Years_since_reg_sc +",
            "Dist_To_Forest_sc +",
            predictor,
            "+",
            interaction_term,
            "+ (1 | Admin_Districts / Subcounty / Village_ID)",
            "+ (1 | Cluster_ID / Group_ID)"
          )
          
          glmmTMB::glmmTMB(
            as.formula(formula_text),
            data = Bushenyi_model_data,
            family = glmmTMB::nbinom2
          )
          
        }
      )
      
    }
  ),
  
  
  tar_target(
    Soroti_H4_trees_models,
    {
      
      purrr::map(
        H4_thresholds,
        function(i) {
          
          predictor <- paste0(
            "NearFar_resid_",
            i,
            "_sc"
          )
          
          interaction_term <- paste0(
            "Exposure_sc:",
            predictor
          )
          
          formula_text <- paste(
            "Trees ~",
            "Exposure_sc +",
            "Years_since_reg_sc +",
            "Exposure_sc:Years_since_reg_sc +",
            "Dist_To_Forest_sc +",
            predictor,
            "+",
            interaction_term,
            "+ (1 | Admin_Districts / Subcounty / Village_ID)",
            "+ (1 | Cluster_ID / Group_ID)"
          )
          
          glmmTMB::glmmTMB(
            as.formula(formula_text),
            data = Soroti_model_data,
            family = glmmTMB::nbinom2
          )
          
        }
      )
      
    }
  ),
  
  
  ################################################################################
  # MODEL CONVERGENCE DIAGNOSTICS
  ################################################################################
  
  
  tar_target(
    mixed_model_convergence,
    bind_rows(
      
      tibble(
        Site = "Bushenyi",
        Outcome = "Density",
        Model = "Base",
        Convergence =
          ifelse(
            is.null(
              Bushenyi_base_density_model@optinfo$conv$lme4$messages
            ),
            "OK",
            "Warning"
          )
      ),
      
      tibble(
        Site = "Soroti",
        Outcome = "Density",
        Model = "Base",
        Convergence =
          ifelse(
            is.null(
              Soroti_base_density_model@optinfo$conv$lme4$messages
            ),
            "OK",
            "Warning"
          )
      ),
      
      tibble(
        Site = "Bushenyi",
        Outcome = "Trees",
        Model = "Base",
        Convergence =
          ifelse(
            is.null(
              Bushenyi_base_trees_model$fit$convergence
            ) ||
              Bushenyi_base_trees_model$fit$convergence == 0,
            "OK",
            "Warning"
          )
      ),
      
      tibble(
        Site = "Soroti",
        Outcome = "Trees",
        Model = "Base",
        Convergence =
          ifelse(
            is.null(
              Soroti_base_trees_model$fit$convergence
            ) ||
              Soroti_base_trees_model$fit$convergence == 0,
            "OK",
            "Warning"
          )
      )
      
    )
  ),
  
  
  ################################################################################
  # BASE VS H4 AIC COMPARISON
  ################################################################################
  
  
  tar_target(
    Bushenyi_base_density_ML,
    update(
      Bushenyi_base_density_model,
      REML = FALSE
    )
  ),
  
  
  tar_target(
    Soroti_base_density_ML,
    update(
      Soroti_base_density_model,
      REML = FALSE
    )
  ),
  
  
  tar_target(
    Bushenyi_base_density_AIC,
    tibble(
      Site = "Bushenyi",
      Outcome = "Density",
      Model = "Base",
      Threshold = NA_integer_,
      AIC = AIC(
        Bushenyi_base_density_ML
      )
    )
  ),
  
  
  tar_target(
    Soroti_base_density_AIC,
    tibble(
      Site = "Soroti",
      Outcome = "Density",
      Model = "Base",
      Threshold = NA_integer_,
      AIC = AIC(
        Soroti_base_density_ML
      )
    )
  ),
  
  
  tar_target(
    Bushenyi_H4_density_AIC,
    tibble(
      Site = "Bushenyi",
      Outcome = "Density",
      Model = "H4",
      Threshold = H4_thresholds,
      AIC = purrr::map_dbl(
        Bushenyi_H4_density_models,
        AIC
      )
    )
  ),
  
  
  tar_target(
    Soroti_H4_density_AIC,
    tibble(
      Site = "Soroti",
      Outcome = "Density",
      Model = "H4",
      Threshold = H4_thresholds,
      AIC = purrr::map_dbl(
        Soroti_H4_density_models,
        AIC
      )
    )
  ),
  
  
  tar_target(
    Bushenyi_base_trees_AIC,
    tibble(
      Site = "Bushenyi",
      Outcome = "Trees",
      Model = "Base",
      Threshold = NA_integer_,
      AIC = AIC(
        Bushenyi_base_trees_model
      )
    )
  ),
  
  
  tar_target(
    Soroti_base_trees_AIC,
    tibble(
      Site = "Soroti",
      Outcome = "Trees",
      Model = "Base",
      Threshold = NA_integer_,
      AIC = AIC(
        Soroti_base_trees_model
      )
    )
  ),
  
  
  tar_target(
    Bushenyi_H4_trees_AIC,
    tibble(
      Site = "Bushenyi",
      Outcome = "Trees",
      Model = "H4",
      Threshold = H4_thresholds,
      AIC = purrr::map_dbl(
        Bushenyi_H4_trees_models,
        AIC
      )
    )
  ),
  
  
  tar_target(
    Soroti_H4_trees_AIC,
    tibble(
      Site = "Soroti",
      Outcome = "Trees",
      Model = "H4",
      Threshold = H4_thresholds,
      AIC = purrr::map_dbl(
        Soroti_H4_trees_models,
        AIC
      )
    )
  ),
  
  
  tar_target(
    mixed_model_AIC_comparison,
    bind_rows(
      
      Bushenyi_base_density_AIC,
      Soroti_base_density_AIC,
      
      Bushenyi_H4_density_AIC,
      Soroti_H4_density_AIC,
      
      Bushenyi_base_trees_AIC,
      Soroti_base_trees_AIC,
      
      Bushenyi_H4_trees_AIC,
      Soroti_H4_trees_AIC
      
    ) %>%
      group_by(
        Site,
        Outcome
      ) %>%
      mutate(
        Delta_AIC =
          AIC - min(AIC)
      ) %>%
      ungroup()
  ),
  
  
  ################################################################################
  # PUBLICATION-READY BASE MODEL RESULTS
  ################################################################################
  
  
  tar_target(
    mixed_model_results,
    bind_rows(
      
      broom.mixed::tidy(
        Bushenyi_base_density_model,
        effects = "fixed"
      ) %>%
        mutate(
          Site = "Bushenyi",
          Outcome = "Density",
          Model = "Base"
        ),
      
      broom.mixed::tidy(
        Soroti_base_density_model,
        effects = "fixed"
      ) %>%
        mutate(
          Site = "Soroti",
          Outcome = "Density",
          Model = "Base"
        ),
      
      broom.mixed::tidy(
        Bushenyi_base_trees_model,
        effects = "fixed"
      ) %>%
        mutate(
          Site = "Bushenyi",
          Outcome = "Trees",
          Model = "Base"
        ),
      
      broom.mixed::tidy(
        Soroti_base_trees_model,
        effects = "fixed"
      ) %>%
        mutate(
          Site = "Soroti",
          Outcome = "Trees",
          Model = "Base"
        )
      
    ) %>%
      select(
        Site,
        Outcome,
        Model,
        everything()
      )
  ),
  
  
  ################################################################################
  # FINAL ERGM TABLES
  ################################################################################
  
  
  tar_target(
    final_ergm_table,
    bind_rows(
      ergm_results
    )
  ),
  
  
  tar_target(
    final_overlap_table,
    bind_rows(
      overlap_results
    )
  ),
  
  
  ################################################################################
  # NATIONAL / SITE SUMMARY STATISTICS + FIGURES 2-5, S1
  #
  # New data target: bushsoroti_raw. Everything else reuses the
  # existing tist_data and districts_sf targets rather than
  # reloading from hardcoded paths.
  ################################################################################
  
  # --------------------------------------------------------------------------
  # New data source
  # --------------------------------------------------------------------------
  
  tar_target(
    bushsoroti_raw_path,
    "Data/TISTDat/BushSoroti_cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered_communicationoptions.csv",
    format = "file"
  ),
  tar_target(bushsoroti_raw, load_bushsoroti_raw(bushsoroti_raw_path)),
  tar_target(bushsoroti_structure_data, prepare_bushsoroti_structure_data(bushsoroti_raw)),
  tar_target(communication_data, prepare_communication_data(bushsoroti_raw)),
  
  
  # --------------------------------------------------------------------------
  # District coverage + national/site summary stats
  # --------------------------------------------------------------------------
  
  tar_target(district_coverage, compute_district_coverage(tist_data, districts_sf)),
  
  tar_target(national_summary_stats, compute_summary_stats(tist_data)),
  tar_target(bushenyi_summary_stats, compute_summary_stats(tist_data, "Bushenyi")),
  tar_target(soroti_summary_stats, compute_summary_stats(tist_data, "Soroti")),
  
  tar_target(national_totals, compute_site_totals(tist_data)),
  tar_target(bushenyi_totals, compute_site_totals(tist_data, "Bushenyi")),
  tar_target(soroti_totals, compute_site_totals(tist_data, "Soroti")),
  
  
  # --------------------------------------------------------------------------
  # Figure 2 -- national 4-panel
  # --------------------------------------------------------------------------
  
  tar_target(groups_per_farmer_national, summarise_groups_per_farmer_national(tist_data)),
  tar_target(farmers_per_group_national, summarise_farmers_per_group_national(tist_data)),
  tar_target(locations_per_farmer_national, summarise_locations_per_farmer_national(tist_data)),
  tar_target(groves_per_village_national, summarise_groves_per_village_national(tist_data)),
  
  tar_target(
    farmer_group_plot,
    plot_count_bar(groups_per_farmer_national, "n_groups_cat", "Number of Groups per Farmer", "Percentage of Farmers")
  ),
  tar_target(
    group_farmer_plot,
    plot_count_bar(farmers_per_group_national, "n_farmers_cat", "Number of Farmers per Group", "Percentage of Groups")
  ),
  tar_target(
    farmer_location_plot,
    plot_count_bar(locations_per_farmer_national, "location_bin", "Number of Grove Locations per Farmer",
                   "Percentage of Farmers", fill = "darkgreen")
  ),
  tar_target(
    location_farmer_plot,
    plot_count_bar(groves_per_village_national, "farmer_bin", "Number of Farmer Groves per Village",
                   "Percentage of Villages", fill = "darkgreen")
  ),
  
  tar_target(
    figure2_national,
    combine_figure2_national(farmer_group_plot, group_farmer_plot, farmer_location_plot, location_farmer_plot)
  ),
  tar_target(
    figure2_national_file,
    save_ggplot(figure2_national, "Output/Manuscript 3 graphs/manuscript_new_plots/Figure2_National_Plot.png",
                width = 183 / 25.4, height = 200 / 25.4, dpi = 600),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Figure 3 -- subcounty penetration + cluster concentration
  # --------------------------------------------------------------------------
  
  tar_target(site_colours, get_site_colours()),
  
  tar_target(farmer_penetration_by_subcounty, compute_farmer_penetration_by_subcounty(bushsoroti_structure_data)),
  tar_target(
    soroti_penetration_plot,
    plot_penetration(farmer_penetration_by_subcounty, "Soroti", threshold = 0.1, site_colours = site_colours)
  ),
  tar_target(
    bushenyi_penetration_plot,
    plot_penetration(farmer_penetration_by_subcounty, "Bushenyi", threshold = 0.5, site_colours = site_colours)
  ),
  
  tar_target(cluster_summary, compute_cluster_summary(bushsoroti_structure_data)),
  tar_target(
    soroti_cluster_plot_data,
    prepare_cluster_plot_data(cluster_summary, "Soroti", threshold = 1, force_other_last = FALSE)
  ),
  tar_target(
    bushenyi_cluster_plot_data,
    prepare_cluster_plot_data(cluster_summary, "Bushenyi", threshold = 2, force_other_last = TRUE)
  ),
  tar_target(soroti_cluster_plot, plot_cluster_bar(soroti_cluster_plot_data, site_colours)),
  tar_target(bushenyi_cluster_plot, plot_cluster_bar(bushenyi_cluster_plot_data, site_colours)),
  
  tar_target(
    figure3_penetration_cluster,
    combine_figure3(soroti_penetration_plot, bushenyi_penetration_plot, soroti_cluster_plot, bushenyi_cluster_plot)
  ),
  tar_target(
    figure3_penetration_cluster_file,
    save_ggplot(figure3_penetration_cluster, "Output/Manuscript 3 graphs/manuscript_new_plots/Figure3_Penetration_Cluster_Plot.png",
                width = 220 / 25.4, height = 220 / 25.4, dpi = 600),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Figure 4 -- site structure comparison
  # --------------------------------------------------------------------------
  
  tar_target(groups_per_farmer_by_site, summarise_groups_per_farmer_by_site(bushsoroti_structure_data)),
  tar_target(farmers_per_group_by_site, summarise_farmers_per_group_by_site(bushsoroti_structure_data)),
  tar_target(locations_per_farmer_by_site, summarise_locations_per_farmer_by_site(bushsoroti_structure_data)),
  tar_target(groves_per_village_by_site, summarise_groves_per_village_by_site(bushsoroti_structure_data)),
  
  tar_target(
    farmer_group_plot_bysite,
    plot_dodged_bar(groups_per_farmer_by_site, "n_groups_cat", "Number of Groups per Farmer",
                    "Percentage of Farmers", site_colours, width = 0.8)
  ),
  tar_target(
    group_farmer_plot_bysite,
    plot_dodged_bar(farmers_per_group_by_site, "n_farmers_cat", "Number of Farmers per Group",
                    "Percentage of Groups", site_colours, width = 0.75)
  ),
  tar_target(
    farmer_location_plot_bysite,
    plot_dodged_bar(locations_per_farmer_by_site, "n_locations_cat", "Number of Groves per Farmer",
                    "Percentage of Farmers", site_colours, width = 0.8)
  ),
  tar_target(
    location_farmer_plot_bysite,
    plot_dodged_bar(groves_per_village_by_site, "n_farmers_bin", "Number of Groves per Village",
                    "Percentage of Villages", site_colours, width = 0.75)
  ),
  
  tar_target(
    figure4_site_structure,
    combine_figure4(farmer_group_plot_bysite, group_farmer_plot_bysite,
                    farmer_location_plot_bysite, location_farmer_plot_bysite)
  ),
  tar_target(
    figure4_site_structure_file,
    save_ggplot(figure4_site_structure, "Output/Manuscript 3 graphs/manuscript_new_plots/Figure4_Site_Structure_Plot.png",
                width = 240 / 25.4, height = 220 / 25.4, dpi = 600),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Figure S1 -- communication channels
  # --------------------------------------------------------------------------
  
  tar_target(comm_vars, get_comm_vars()),
  tar_target(subcounty_comm_data, prepare_subcounty_comm_data(communication_data, comm_vars)),
  tar_target(comm_channel_summary, summarise_comm_channels(subcounty_comm_data, comm_vars)),
  tar_target(figureS1_communication, plot_comm_channels(comm_channel_summary)),
  tar_target(
    figureS1_communication_file,
    save_ggplot(figureS1_communication, "Output/Manuscript 3 graphs/manuscript_new_plots/FigureS1_Communication_Plot.png",
                width = 183 / 25.4, height = 120 / 25.4, dpi = 600),
    format = "file"
  ),
  
  
  # --------------------------------------------------------------------------
  # Figure 5 -- ridgeline plots (tree count, area, density)
  # --------------------------------------------------------------------------
  
  tar_target(tist_data_with_site, add_site_column(tist_data)),
  
  tar_target(
    trees_plot_data,
    tist_data_with_site |> dplyr::filter(!is.na(Trees), Trees > 0)
  ),
  tar_target(trees_ridgeline_stats, compute_ridgeline_stats(trees_plot_data, "Trees")),
  tar_target(
    ridgeline_tree_plot,
    plot_ridgeline(trees_plot_data, "Trees", trees_ridgeline_stats,
                   "Number of trees per farmer (log scale)",
                   digits_med = NULL, digits_mean = 1, digits_p90 = 0)
  ),
  
  tar_target(
    area_plot_data,
    tist_data_with_site |> dplyr::filter(!is.na(Area_Ha), Area_Ha > 0)
  ),
  tar_target(min_area, min(area_plot_data$Area_Ha)),
  tar_target(area_ridgeline_stats, compute_ridgeline_stats(area_plot_data, "Area_Ha")),
  tar_target(
    ridgeline_area_plot,
    plot_ridgeline(area_plot_data, "Area_Ha", area_ridgeline_stats,
                   "Farm area (ha, log scale)",
                   digits_med = 2, digits_mean = 2, digits_p90 = 2, x_limits = c(min_area, NA))
  ),
  
  tar_target(
    density_plot_data,
    tist_data_with_site |> dplyr::filter(!is.na(Density_winsor99), Density_winsor99 > 0)
  ),
  tar_target(min_density, min(density_plot_data$Density_winsor99, na.rm = TRUE)),
  tar_target(density_ridgeline_stats, compute_ridgeline_stats(density_plot_data, "Density_winsor99")),
  tar_target(
    ridgeline_density_plot,
    plot_ridgeline(density_plot_data, "Density_winsor99", density_ridgeline_stats,
                   "Planted tree density, 99th percentile capped (trees per hectare, log scale)",
                   digits_med = 0, digits_mean = 0, digits_p90 = 0, trim = TRUE, from = min_density)
  ),
  
  tar_target(
    figure5_ridgeline,
    combine_ridgeline_figure(ridgeline_tree_plot, ridgeline_area_plot, ridgeline_density_plot)
  ),
  tar_target(
    figure5_ridgeline_file,
    save_ggplot(figure5_ridgeline, "Output/Manuscript 3 graphs/manuscript_new_plots/Figure5_Ridgeline_Plot.png",
                width = 180 / 25.4, height = 260 / 25.4, dpi = 600),
    format = "file"
  )
  
)
