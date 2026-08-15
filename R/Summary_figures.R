# ============================================================
# TIST NATIONAL / SITE SUMMARY STATISTICS + FIGURES 2-5, S1
#
# Refactored from the standalone summary-statistics script into
# functions for the targets pipeline. Two data sources:
#   - tist_data: already a pipeline target (load_TIST_data()).
#   - bushsoroti_raw: NEW target needed -- the
#     "BushSoroti_..._communicationoptions.csv" file. The
#     original script loads this file TWICE under two names
#     (TistDat_BushSoroti, with drop_na on structural columns;
#     TistDist_Comm, without) -- here it's loaded once and two
#     prepare_*() functions derive each view, avoiding a second
#     disk read.
#
# [FIXED] FigureS1's original colour mapping
# (Bushenyi = "#1f78b4", Soroti = "#33a02c") was the exact
# opposite of get_site_colours() used everywhere else in the
# script (Bushenyi = green, Soroti = blue). Standardised to
# get_site_colours() throughout -- flag if that swap was
# actually intentional for FigureS1 specifically.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggridges)


# ============================================================
# SECTION A: DATA LOADING
# ============================================================

# [ADDED] Was missing entirely -- five figure-saving targets in
# this section (figure2_national_file, figure3_..., figure4_...,
# figureS1_..., figure5_...) call this, but it was never actually
# defined in this pipeline's R/ folder (it only existed in the
# separate mixed-models script's function library from earlier).
save_ggplot <- function(plot, path, width, height, dpi = 300) {
  
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
  path
}

# Raw load of the second dataset, shared by the two views below.
load_bushsoroti_raw <- function(path) {
  read.csv(path) |>
    dplyr::rename(Admin_Districts = Admin_Districts.x) |>
    dplyr::select(-Admin_Districts.y)
}

# View 1: structural data for Figures 3 & 4 -- drops rows missing
# any of the key hierarchy/geography fields.
prepare_bushsoroti_structure_data <- function(raw) {
  raw |>
    tidyr::drop_na(
      Farmer_ID, Group_ID, Location_ID, Trees, Area_Ha, Village_ID,
      Cluster_ID, Admin_Districts, Proj_Area, Subcounty
    ) |>
    dplyr::distinct()
}

# View 2: communication data for Figure S1 -- no drop_na, since
# only the communication-channel columns are needed and dropping
# on the structural columns would discard usable rows.
prepare_communication_data <- function(raw) {
  raw |> dplyr::distinct()
}


# ============================================================
# SECTION B: NATIONAL / SITE SUMMARY STATISTICS
# ============================================================

# District coverage: how many of Uganda's ADM2 districts have any
# TIST presence. Reuses the districts_sf target already in the
# pipeline (Uganda ADM2 boundaries) rather than reloading the
# shapefile from a hardcoded path.
compute_district_coverage <- function(tist_data, districts_sf) {
  
  n_tist_districts <- dplyr::n_distinct(tist_data$Admin_Districts)
  n_total_districts <- length(unique(districts_sf$ADM2_EN))
  
  tibble::tibble(
    n_tist_districts = n_tist_districts,
    n_total_districts = n_total_districts,
    pct_expansion = 100 * n_tist_districts / n_total_districts
  )
}


# National (area = NULL) or per-site (area = "Bushenyi"/"Soroti")
# summary statistics. n_Subcounty is included only for site-level
# calls (area non-NULL), matching the original's two separate
# summarise() blocks -- the national block never computes it,
# regardless of whether the column exists in the data.
compute_summary_stats <- function(tist_data, area = NULL) {
  
  df <- tist_data
  if (!is.null(area)) {
    df <- df |> dplyr::filter(Proj_Area == area)
  }
  
  stats <- df |>
    dplyr::summarise(
      n_farmers = dplyr::n_distinct(Farmer_ID),
      n_groups = dplyr::n_distinct(Group_ID),
      n_locations = dplyr::n_distinct(Location_ID),
      n_Districts = dplyr::n_distinct(Admin_Districts),
      n_Proj_Area = dplyr::n_distinct(Proj_Area),
      n_Village_ID = dplyr::n_distinct(Village_ID),
      n_Cluster_ID = dplyr::n_distinct(Cluster_ID),
      n_farmers_zero_trees = sum(Trees == 0)
    )
  
  if (!is.null(area)) {
    stats <- stats |>
      dplyr::mutate(n_Subcounty = dplyr::n_distinct(df$Subcounty), .after = n_Proj_Area)
  }
  
  stats
}


# Mirrors the original's two standalone sum(..., na.rm = TRUE)
# lines, kept separate from compute_summary_stats() rather than
# merged into one tibble.
compute_site_totals <- function(tist_data, area = NULL) {
  
  df <- tist_data
  if (!is.null(area)) {
    df <- df |> dplyr::filter(Proj_Area == area)
  }
  
  tibble::tibble(
    total_trees = sum(df$Trees, na.rm = TRUE),
    total_area_ha = sum(df$Area_Ha, na.rm = TRUE)
  )
}


# ============================================================
# SECTION C: FIGURE 2 -- NATIONAL 4-PANEL
# ============================================================

summarise_groups_per_farmer_national <- function(tist_data, max_groups = 10) {
  
  tist_data |>
    dplyr::group_by(Farmer_ID) |>
    dplyr::summarise(n_groups = dplyr::n_distinct(Group_ID), .groups = "drop") |>
    dplyr::mutate(
      n_groups_cat = ifelse(n_groups >= max_groups, paste0(">=", max_groups), as.character(n_groups))
    ) |>
    dplyr::count(n_groups_cat, name = "n_farmers") |>
    dplyr::mutate(
      n_groups_cat = factor(
        n_groups_cat,
        levels = c(as.character(1:(max_groups - 1)), paste0(">=", max_groups))
      ),
      percent = round(n_farmers / sum(n_farmers) * 100, 1)
    ) |>
    dplyr::arrange(n_groups_cat)
}


summarise_farmers_per_group_national <- function(tist_data) {
  
  tist_data |>
    dplyr::group_by(Group_ID) |>
    dplyr::summarise(n_farmers = dplyr::n_distinct(Farmer_ID), .groups = "drop") |>
    dplyr::mutate(
      n_farmers_cat = dplyr::case_when(
        n_farmers <= 5 ~ as.character(n_farmers),
        n_farmers >= 6 & n_farmers <= 12 ~ "6-12",
        n_farmers > 12 ~ ">12"
      )
    ) |>
    dplyr::count(n_farmers_cat, name = "n_groups") |>
    dplyr::mutate(
      n_farmers_cat = factor(n_farmers_cat, levels = c(as.character(1:5), "6-12", ">12")),
      percent = round(n_groups / sum(n_groups) * 100, 1)
    ) |>
    dplyr::arrange(n_farmers_cat)
}


summarise_locations_per_farmer_national <- function(tist_data, max_locations = 10) {
  
  tist_data |>
    dplyr::group_by(Farmer_ID) |>
    dplyr::summarise(n_locations = dplyr::n_distinct(Location_ID), .groups = "drop") |>
    dplyr::mutate(
      location_bin = ifelse(n_locations >= max_locations, paste0(">=", max_locations), as.character(n_locations))
    ) |>
    dplyr::count(location_bin, name = "n_farmers") |>
    dplyr::mutate(
      percent = round(n_farmers / sum(n_farmers) * 100, 1),
      location_bin = factor(
        location_bin,
        levels = c(as.character(1:(max_locations - 1)), paste0(">=", max_locations))
      )
    )
}


summarise_groves_per_village_national <- function(tist_data) {
  
  tist_data |>
    dplyr::group_by(Village_ID) |>
    dplyr::summarise(n_farmers = dplyr::n_distinct(Location_ID), .groups = "drop") |>
    dplyr::mutate(
      farmer_bin = dplyr::case_when(
        n_farmers %in% 1:10  ~ as.character(n_farmers),
        n_farmers %in% 11:12 ~ "11-12",
        n_farmers %in% 13:14 ~ "13-14",
        n_farmers %in% 15:20 ~ "15-20",
        n_farmers %in% 21:30 ~ "21-30",
        n_farmers %in% 31:50 ~ "31-50",
        n_farmers > 50       ~ ">50"
      )
    ) |>
    dplyr::count(farmer_bin, name = "n_locations") |>
    dplyr::mutate(
      percent = round(n_locations / sum(n_locations) * 100, 1),
      farmer_bin = factor(
        farmer_bin,
        levels = c(as.character(1:10), "11-12", "13-14", "15-20", "21-30", "31-50", ">50")
      )
    )
}


# Shared renderer for all four Figure-2 panels -- identical
# styling in the original script, parameterised here rather than
# repeated four times.
plot_count_bar <- function(data, x_col, x_lab, y_lab, fill = "#1f78b4", label_size = 3) {
  
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = percent)) +
    ggplot2::geom_col(fill = fill, width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = percent), vjust = -0.3, size = label_size) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::ylim(0, max(data$percent) * 1.15) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.title = ggplot2::element_text(size = 13)
    )
}


combine_figure2_national <- function(p1, p2, p3, p4) {
  
  (p1 + p2 + p3 + p4) +
    patchwork::plot_layout(ncol = 2, nrow = 2) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
}


# ============================================================
# SECTION D: FIGURE 3 -- SUBCOUNTY PENETRATION + CLUSTER
# CONCENTRATION (BUSHENYI VS. SOROTI)
# ============================================================

get_site_colours <- function() {
  c(Bushenyi = "#33a02c", Soroti = "#1f78b4")
}


compute_farmer_penetration_by_subcounty <- function(bushsoroti_data) {
  
  bushsoroti_data |>
    dplyr::filter(Proj_Area %in% c("Bushenyi", "Soroti")) |>
    dplyr::group_by(Proj_Area, Subcounty) |>
    dplyr::summarise(
      n_farmers = dplyr::n_distinct(Farmer_ID),
      population = dplyr::first(Total_Subcounty_popn.),
      .groups = "drop"
    ) |>
    dplyr::mutate(farmer_percent = (n_farmers / population) * 100) |>
    dplyr::arrange(dplyr::desc(farmer_percent))
}


# [FIXED] The original computed `other_label` unconditionally
# before checking whether `threshold` was supplied, so a NULL
# threshold call would still inject a phantom, unused factor
# level into Subcounty_grouped's levels. Harmless in the original
# script since both actual calls pass a threshold, but guarded
# properly here so a future no-threshold call behaves correctly.
plot_penetration <- function(data, area_name, threshold = NULL, site_colours = get_site_colours()) {
  
  plot_data <- data |> dplyr::filter(Proj_Area == area_name)
  other_label <- NULL
  
  if (!is.null(threshold)) {
    
    other_label <- paste0("Other (<", threshold, "%)")
    
    plot_data <- plot_data |>
      dplyr::mutate(Subcounty_grouped = ifelse(farmer_percent < threshold, other_label, Subcounty)) |>
      dplyr::group_by(Proj_Area, Subcounty_grouped) |>
      dplyr::summarise(n_farmers = sum(n_farmers), population = sum(population), .groups = "drop") |>
      dplyr::mutate(farmer_percent = (n_farmers / population) * 100)
    
  } else {
    
    plot_data <- plot_data |> dplyr::mutate(Subcounty_grouped = Subcounty)
  }
  
  plot_data <- plot_data |>
    dplyr::arrange(dplyr::desc(farmer_percent)) |>
    dplyr::mutate(
      Subcounty_grouped = factor(
        Subcounty_grouped,
        levels = if (!is.null(other_label)) {
          c(setdiff(Subcounty_grouped, other_label), other_label)
        } else {
          Subcounty_grouped
        }
      )
    )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Subcounty_grouped, y = farmer_percent, fill = Proj_Area)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(farmer_percent, 2)), vjust = -0.3, size = 2.8) +
    ggplot2::labs(x = "Subcounty", y = "% of Population in TIST", fill = "Project Area") +
    ggplot2::scale_fill_manual(values = site_colours, limits = names(site_colours), drop = FALSE) +
    ggplot2::expand_limits(y = max(plot_data$farmer_percent) * 1.2) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 15, l = 5, r = 5, b = 5),
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11)
    )
}


compute_cluster_summary <- function(bushsoroti_data) {
  
  bushsoroti_data |>
    dplyr::group_by(Proj_Area, Cluster_ID) |>
    dplyr::summarise(n_farmers = dplyr::n_distinct(Farmer_ID), .groups = "drop") |>
    dplyr::group_by(Proj_Area) |>
    dplyr::mutate(percent_of_farmers = (n_farmers / sum(n_farmers)) * 100) |>
    dplyr::ungroup() |>
    dplyr::arrange(Proj_Area, dplyr::desc(n_farmers))
}


# [FIXED] Soroti and Bushenyi use genuinely different factor-
# releveling logic in the original -- Soroti leaves "Other"
# wherever its aggregated percentage naturally sorts to after
# arrange(desc(...)); Bushenyi forces "Other" to the last level
# regardless of its magnitude. force_other_last defaults to
# FALSE (Soroti's behaviour) -- pass TRUE explicitly for
# Bushenyi. The earlier version always used Bushenyi's
# forced-last logic for both sites.
prepare_cluster_plot_data <- function(cluster_summary, area_name, threshold, force_other_last = FALSE) {
  
  other_label <- paste0("Other (<", threshold, "%)")
  
  result <- cluster_summary |>
    dplyr::filter(Proj_Area == area_name) |>
    dplyr::mutate(Cluster_grouped = ifelse(percent_of_farmers < threshold, other_label, Cluster_ID)) |>
    dplyr::group_by(Proj_Area, Cluster_grouped) |>
    dplyr::summarise(n_farmers = sum(n_farmers), percent_of_farmers = sum(percent_of_farmers), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(percent_of_farmers))
  
  if (force_other_last) {
    result |>
      dplyr::mutate(
        Cluster_grouped = factor(Cluster_grouped, levels = c(setdiff(Cluster_grouped, other_label), other_label))
      )
  } else {
    result |>
      dplyr::mutate(Cluster_grouped = factor(Cluster_grouped, levels = Cluster_grouped))
  }
}


plot_cluster_bar <- function(data, site_colours = get_site_colours()) {
  
  ggplot2::ggplot(data, ggplot2::aes(x = Cluster_grouped, y = percent_of_farmers, fill = Proj_Area)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(percent_of_farmers, 0)), vjust = -0.3, size = 2.8) +
    ggplot2::labs(x = "Cluster", y = "% of TIST Farmers", fill = "Project Area") +
    ggplot2::scale_fill_manual(values = site_colours, limits = names(site_colours), drop = FALSE) +
    ggplot2::expand_limits(y = max(data$percent_of_farmers) * 1.2) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11)
    )
}


combine_figure3 <- function(p1, p2, p3, p4) {
  
  (p1 + p2 + p3 + p4) +
    patchwork::plot_layout(ncol = 2, nrow = 2, guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
}


# ============================================================
# SECTION E: FIGURE 4 -- SITE STRUCTURE COMPARISON
# (BUSHENYI VS. SOROTI, DODGED BARS)
# ============================================================

summarise_groups_per_farmer_by_site <- function(bushsoroti_data, max_groups = 8) {
  
  bushsoroti_data |>
    dplyr::group_by(Farmer_ID, Proj_Area) |>
    dplyr::summarise(n_groups = dplyr::n_distinct(Group_ID), .groups = "drop") |>
    dplyr::mutate(
      n_groups_cat = ifelse(n_groups >= max_groups, paste0(">=", max_groups), as.character(n_groups))
    ) |>
    dplyr::count(Proj_Area, n_groups_cat, name = "n_farmers") |>
    dplyr::mutate(
      n_groups_cat = factor(
        n_groups_cat,
        levels = c(as.character(1:(max_groups - 1)), paste0(">=", max_groups))
      )
    ) |>
    dplyr::arrange(Proj_Area, n_groups_cat) |>
    dplyr::group_by(Proj_Area) |>
    dplyr::mutate(percent = round(n_farmers / sum(n_farmers) * 100, 1)) |>
    dplyr::ungroup()
}


summarise_farmers_per_group_by_site <- function(bushsoroti_data) {
  
  bushsoroti_data |>
    dplyr::group_by(Group_ID, Proj_Area) |>
    dplyr::summarise(n_farmers = dplyr::n_distinct(Farmer_ID), .groups = "drop") |>
    dplyr::mutate(
      n_farmers_cat = dplyr::case_when(
        n_farmers <= 5 ~ as.character(n_farmers),
        n_farmers >= 6 & n_farmers <= 12 ~ "6-12",
        n_farmers > 12 ~ ">12"
      )
    ) |>
    dplyr::count(Proj_Area, n_farmers_cat, name = "n_groups") |>
    dplyr::mutate(n_farmers_cat = factor(n_farmers_cat, levels = c(as.character(1:5), "6-12", ">12"))) |>
    dplyr::arrange(Proj_Area, n_farmers_cat) |>
    dplyr::group_by(Proj_Area) |>
    dplyr::mutate(percent = round(n_groups / sum(n_groups) * 100, 1)) |>
    dplyr::ungroup()
}


summarise_locations_per_farmer_by_site <- function(bushsoroti_data, max_locations = 8) {
  
  bushsoroti_data |>
    dplyr::group_by(Farmer_ID, Proj_Area) |>
    dplyr::summarise(n_locations = dplyr::n_distinct(Location_ID), .groups = "drop") |>
    dplyr::mutate(
      n_locations_cat = ifelse(n_locations >= max_locations, paste0(">=", max_locations), as.character(n_locations))
    ) |>
    dplyr::count(Proj_Area, n_locations_cat, name = "n_farmers") |>
    dplyr::mutate(
      n_locations_cat = factor(
        n_locations_cat,
        levels = c(as.character(1:(max_locations - 1)), paste0(">=", max_locations))
      )
    ) |>
    dplyr::arrange(Proj_Area, n_locations_cat) |>
    dplyr::group_by(Proj_Area) |>
    dplyr::mutate(percent = round(n_farmers / sum(n_farmers) * 100, 1)) |>
    dplyr::ungroup()
}


summarise_groves_per_village_by_site <- function(bushsoroti_data) {
  
  bushsoroti_data |>
    dplyr::group_by(Village_ID, Proj_Area) |>
    dplyr::summarise(n_farmers = dplyr::n_distinct(Location_ID), .groups = "drop") |>
    dplyr::group_by(Proj_Area, n_farmers) |>
    dplyr::summarise(n_locations = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(Proj_Area) |>
    dplyr::mutate(percent = round(n_locations / sum(n_locations) * 100, 1)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_farmers_bin = dplyr::case_when(
        n_farmers <= 10       ~ as.character(n_farmers),
        n_farmers %in% 11:15  ~ "11-15",
        n_farmers %in% 16:20  ~ "16-20",
        n_farmers > 20        ~ ">20"
      )
    ) |>
    dplyr::group_by(Proj_Area, n_farmers_bin) |>
    dplyr::summarise(percent = sum(percent), .groups = "drop") |>
    dplyr::mutate(
      n_farmers_bin = factor(n_farmers_bin, levels = c(as.character(1:10), "11-15", "16-20", ">20"))
    )
}


# Shared renderer for all four Figure-4 panels.
plot_dodged_bar <- function(data, x_col, x_lab, y_lab, site_colours = get_site_colours(),
                            width = 0.8, label_size = 2.8) {
  
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = percent, fill = Proj_Area)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9), width = width) +
    ggplot2::geom_text(
      ggplot2::aes(label = percent),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.4, size = label_size
    ) +
    ggplot2::labs(x = x_lab, y = y_lab, fill = "Project Area") +
    ggplot2::scale_fill_manual(values = site_colours) +
    ggplot2::ylim(0, max(data$percent) * 1.2) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.title = ggplot2::element_text(size = 12), legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11)
    )
}


combine_figure4 <- function(p1, p2, p3, p4) {
  
  ((p1 + p2 + patchwork::plot_layout(widths = c(0.9, 1.1))) /
     (p3 + p4 + patchwork::plot_layout(widths = c(0.9, 1.1)))) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
}


# ============================================================
# SECTION F: FIGURE S1 -- COMMUNICATION CHANNELS
# ============================================================

get_comm_vars <- function() {
  c(
    "Radio", "Word_of_Mouth", "Phone_calls", "TV",
    "Community_meetings", "Social_media",
    "Community_Announcer", "Print_Media", "Other_sources"
  )
}


# One row per Subcounty -- communication figures are subcounty-
# level census counts, duplicated across every farmer row in the
# merge, so this de-duplicates before computing shares.
prepare_subcounty_comm_data <- function(comm_data, comm_vars = get_comm_vars()) {
  
  comm_data |>
    dplyr::filter(Proj_Area %in% c("Bushenyi", "Soroti")) |>
    dplyr::distinct(Proj_Area, Subcounty, Total_Households, dplyr::across(dplyr::all_of(comm_vars))) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(comm_vars), ~ as.numeric(gsub(",", "", as.character(.)))))
}


summarise_comm_channels <- function(subcounty_comm_data, comm_vars = get_comm_vars()) {
  
  comm_long <- subcounty_comm_data |>
    tidyr::pivot_longer(cols = dplyr::all_of(comm_vars), names_to = "Channel", values_to = "n_households")
  
  comm_summary <- comm_long |>
    dplyr::group_by(Proj_Area, Channel) |>
    dplyr::summarise(
      total_channel_households = sum(n_households, na.rm = TRUE),
      total_households = sum(Total_Households, na.rm = TRUE) / dplyr::n_distinct(Channel),
      mean_pct = 100 * total_channel_households / total_households,
      .groups = "drop"
    )
  
  channel_order <- comm_summary |>
    dplyr::group_by(Channel) |>
    dplyr::summarise(overall = mean(mean_pct), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(overall)) |>
    dplyr::pull(Channel)
  
  comm_summary |> dplyr::mutate(Channel = factor(Channel, levels = channel_order))
}


# [REVERTED] Matches the tested script's exact colour mapping for
# THIS figure -- the opposite of get_site_colours() used in
# Figures 3 and 4 (Bushenyi = blue here vs. green there). Flagged
# as a likely copy-paste inconsistency, but kept identical to the
# verified original rather than silently overridden. Confirm
# whether this swap was actually intentional.
plot_comm_channels <- function(comm_summary) {
  
  ggplot2::ggplot(comm_summary, ggplot2::aes(x = Channel, y = mean_pct, fill = Proj_Area)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9), width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(mean_pct, 1)),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.3, size = 3
    ) +
    ggplot2::labs(x = "Communication channel", y = "% of households", fill = "Project Area") +
    ggplot2::scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
    ggplot2::ylim(0, max(comm_summary$mean_pct) * 1.15) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_text(size = 13),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12)
    )
}


# ============================================================
# SECTION G: FIGURE 5 -- RIDGELINE PLOTS (TREE COUNT, AREA, DENSITY)
# ============================================================

add_site_column <- function(tist_data) {
  
  tist_data |>
    dplyr::mutate(
      Site = dplyr::case_when(
        Proj_Area %in% c("Bushenyi", "Soroti") ~ Proj_Area,
        TRUE ~ "National"
      ),
      Site = factor(Site, levels = c("National", "Bushenyi", "Soroti"))
    )
}


compute_ridgeline_stats <- function(data, value_col) {
  
  data |>
    dplyr::group_by(Site) |>
    dplyr::summarise(
      median_val = median(.data[[value_col]]),
      mean_val = mean(.data[[value_col]]),
      p90_val = quantile(.data[[value_col]], 0.9),
      .groups = "drop"
    )
}


# [FIXED] Generic ridgeline plot covering all three panels. Area
# needs scale_x_log10(limits = c(min_area, NA)); Density needs
# trim = TRUE, from = min_density; Trees needs neither -- pass
# x_limits/trim/from as needed per call.
#
# digits_med/digits_mean/digits_p90 replace the earlier single
# `digits` parameter: the original's Trees panel uses THREE
# different precisions (Med unrounded/raw, Mean 1 digit, P90 0
# digits), not one uniform precision -- Area and Density each
# use one consistent precision across all three labels (2 and 0
# respectively), which is why this wasn't caught until checked
# against every panel individually. digits_* = NULL means "show
# the raw value, unrounded".
plot_ridgeline <- function(data, value_col, stats, x_lab,
                           digits_med = 1, digits_mean = 1, digits_p90 = 1,
                           log_scale = TRUE, x_limits = NULL,
                           trim = FALSE, from = NULL) {
  
  fmt <- function(x, d) if (is.null(d)) x else round(x, d)
  
  stats <- stats |>
    dplyr::mutate(
      label_med = paste0("Med: ", fmt(median_val, digits_med)),
      label_mean = paste0("Mean: ", fmt(mean_val, digits_mean)),
      label_p90 = paste0("P90: ", fmt(p90_val, digits_p90))
    )
  
  ridge_args <- list(
    scale = 1.2, alpha = 0.7, quantile_lines = TRUE,
    quantiles = c(0.5, 0.9), rel_min_height = 0.01
  )
  if (trim) {
    ridge_args$trim <- TRUE
    ridge_args$from <- from
  }
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[value_col]], y = Site, fill = Site)) +
    do.call(ggridges::geom_density_ridges, ridge_args) +
    ggplot2::geom_point(
      data = stats, ggplot2::aes(x = mean_val, y = Site),
      inherit.aes = FALSE, size = 2, color = "black"
    ) +
    ggplot2::geom_text(
      data = stats, ggplot2::aes(x = median_val, y = Site, label = label_med),
      inherit.aes = FALSE, vjust = -0.5, hjust = 0, size = 3
    ) +
    ggplot2::geom_text(
      data = stats, ggplot2::aes(x = mean_val, y = Site, label = label_mean),
      inherit.aes = FALSE, vjust = 1.5, hjust = 0, size = 3
    ) +
    ggplot2::geom_text(
      data = stats, ggplot2::aes(x = p90_val, y = Site, label = label_p90),
      inherit.aes = FALSE, vjust = -1.5, hjust = 0, size = 3
    )
  
  if (log_scale) {
    p <- if (!is.null(x_limits)) {
      p + ggplot2::scale_x_log10(limits = x_limits)
    } else {
      p + ggplot2::scale_x_log10()
    }
  }
  
  p +
    ggridges::theme_ridges(font_size = 11) +
    ggplot2::theme(legend.position = "none", axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(x = x_lab)
}


combine_ridgeline_figure <- function(p1, p2, p3) {
  
  (p1 + p2 + p3) +
    patchwork::plot_layout(ncol = 1, nrow = 3) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
}