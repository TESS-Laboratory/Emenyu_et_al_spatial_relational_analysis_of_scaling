# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(magrittr)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "ergm.multi", "network",
               "geosphere", "ergm", "sna", "statnet",
               "patchwork", "modelsummary", "terra","sf", "scales" )
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
      "Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv"
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
  
  # ESA WorldCover 2020 tiles covering Uganda
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
  
  # Persistent land-cover mosaic on disk
  tar_target(
    landcover_raster_file,
    build_landcover_mosaic_file(
      files = landcover_files,
      boundary = Uganda,
      output = "Data/landcover/ESA_WorldCover_2020/Uganda_WorldCover_2020_mosaic.tif"
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
  
  # ============================================================
  # ELIGIBLE LAND / FOREST EXCLUSION
  # ============================================================
  
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
  
  # ------------------------------------------------------------
  # Forest reserves prepared as a persistent vector file
  # in the land-cover raster CRS for rasterisation
  # ------------------------------------------------------------
  
  tar_target(
    forest_reserves_vector_file,
    prepare_forest_reserves_for_exclusion_file(
      forest_reserves,
      landcover_raster_file
    ),
    format = "file"
  ),
  
  # ------------------------------------------------------------
  # Forest reserves prepared as a persistent vector file
  # in UTM Zone 36N for distance calculations
  # ------------------------------------------------------------
  
  tar_target(
    forest_reserves_projected_file,
    prepare_forest_reserves_file(
      forest_reserves,
      target_crs = 32636
    ),
    format = "file"
  ),
  
  # ------------------------------------------------------------
  # Rasterise forest reserves
  # ------------------------------------------------------------
  
  tar_target(
    forest_mask_file,
    create_forest_exclusion_mask_file(
      forest_reserves_vector_file,
      landcover_raster_file
    ),
    format = "file"
  ),
  
  # ------------------------------------------------------------
  # Remove forest reserves from eligible land
  # ------------------------------------------------------------
  
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
      "Buhweju", "Bushenyi", "Ibanda", "Kamwenge", "Kasese",
      "Kitagwenda", "Kyenjojo", "Mitooma", "Ntungamo",
      "Rubirizi", "Rwampara", "Sheema"
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Soroti district names
  # --------------------------------------------------------------------------
  
  tar_target(
    Soroti_names,
    c(
      "Alebtong", "Amuria", "Kalaki", "Kapelebyong", "Serere", "Soroti"
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
  # [FIXED] Was mask_project_area(eligible_area_file, Bushenyi_area),
  # storing a live SpatRaster as the target's value. SpatRaster
  # holds an external C++ pointer that goes invalid the moment
  # targets serialises it to disk and reloads it -- this is what
  # caused "external pointer is not valid" in random_Soroti_points
  # downstream. Switched to the file-persisting wrapper, matching
  # every other raster target in this pipeline.
  # prepare_random_ag_points() already accepts a file path
  # directly, so nothing downstream needs to change.
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
  # [FIXED] Same live-SpatRaster -> file-path correction.
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
  # [BUG FIX] Removed the stray "tar" that was glued onto the end
  # of the print(summary(...)) line -- this broke the parse of
  # this target's command and is almost certainly what caused the
  # "argument 45 is empty" error.
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
  # ERGM OVERLAP SIMULATIONS -- DISABLED
  #
  # nsim = 1000 per site x distance-threshold combination is
  # expensive and results already exist outside this pipeline.
  # Commented out along with simulation_diagnostics,
  # overlap_results, and final_overlap_table below (all depend on
  # this), since leaving any of those active while this is
  # disabled would error with a missing-target reference.
  ################################################################################
  
  # tar_target(
  #   overlap_simulations,
  #   evaluate_overlap(
  #     site_models,
  #     nsim = 1000
  #   ),
  #   pattern = map(site_models)
  # ),
  
  
  ################################################################################
  # SIMULATION DIAGNOSTICS -- DISABLED (depends on overlap_simulations)
  ################################################################################
  
  # tar_target(
  #   simulation_diagnostics,
  #   {
  #
  #     cat(
  #       "\n==============================\n",
  #       "Site:",
  #       overlap_simulations$site,
  #       "\nDistance:",
  #       overlap_simulations$distance_threshold,
  #       "m\n",
  #       "==============================\n"
  #     )
  #
  #     print(
  #       overlap_simulations$observed
  #     )
  #
  #     print(
  #       overlap_simulations$p.values
  #     )
  #
  #     overlap_simulations
  #
  #   },
  #   pattern = map(overlap_simulations)
  # ),
  
  
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
  # PUBLICATION-READY OVERLAP RESULTS -- DISABLED (depends on overlap_simulations)
  ################################################################################
  
  # tar_target(
  #   overlap_results,
  #   {
  #
  #     tibble(
  #
  #       Site =
  #         overlap_simulations$site,
  #
  #       Distance =
  #         overlap_simulations$distance_threshold,
  #
  #       Statistic =
  #         c(
  #           "gcor",
  #           "Overlap"
  #         ),
  #
  #       Observed =
  #         c(
  #           overlap_simulations$observed$gcor,
  #           overlap_simulations$observed$overlap
  #         ),
  #
  #       Simulated_mean =
  #         c(
  #           mean(
  #             overlap_simulations$simulated$gcor
  #           ),
  #           mean(
  #             overlap_simulations$simulated$overlap
  #           )
  #         ),
  #
  #       Simulated_SD =
  #         c(
  #           sd(
  #             overlap_simulations$simulated$gcor
  #           ),
  #           sd(
  #             overlap_simulations$simulated$overlap
  #           )
  #         ),
  #
  #       P_value =
  #         c(
  #           overlap_simulations$p.values$gcor,
  #           overlap_simulations$p.values$overlap
  #         )
  #     )
  #
  #   },
  #   pattern = map(overlap_simulations)
  # ),
  
  
  ################################################################################
  # FINAL ERGM TABLES
  ################################################################################
  
  tar_target(
    final_ergm_table,
    bind_rows(
      ergm_results
    )
  )
  
  # [DISABLED -- depends on overlap_results, which is disabled above]
  # ,
  # tar_target(
  #   final_overlap_table,
  #   bind_rows(
  #     overlap_results
  #   )
  # )
  
)
