################################################################################
# TIST Uganda — Descriptive Statistics & Manuscript Figures
# (Bushenyi & Soroti: farmer/group/location distributions, district
#  penetration, cluster composition, communication channels, ridgelines)
################################################################################

## ---- Packages ------------------------------------------------------------
# install.packages("ggridges")  # run once
library(ggridges)
library(tidyverse)
library(patchwork)
library(sf)
library(flextable)
library(officer)


################################################################################
# 1. LOAD DATA
################################################################################

TistDist_Dat <- read.csv("Data/TISTDat/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv") %>%
  drop_na() %>%
  distinct()

head(TistDist_Dat)
nrow(TistDist_Dat)
n_distinct(TistDist_Dat$Farmer_ID) 

TistDist_Bush <- TistDist_Dat %>%
  filter(Proj_Area %in% c("Bushenyi"))
head(TistDist_Bush)

TistDist_Soroti <- TistDist_Dat %>%
  filter(Proj_Area %in% c("Soroti"))
head(TistDist_Soroti)

Tist_Bush_Soroti <- TistDist_Dat %>%
  filter(Proj_Area %in% c("Soroti", "Bushenyi"))
head(Tist_Bush_Soroti)


################################################################################
# 2. NATIONAL DISTRICT COVERAGE
################################################################################

Admin_Districts <- st_read("C:/Users/ae474/OneDrive - University of Exeter/PhD work in progress/Data/R Large Data/UGshapefiles/uga_admbnda_adm2_ubos_20200824.shp")
head(Admin_Districts)

# Count unique district names
length(unique(Admin_Districts$ADM2_EN))    ## 135
length(unique(TistDist_Dat$Admin_Districts)) ## 57

## Proportion of districts with TIST presence
Percentate_expansion_by_district <- length(unique(TistDist_Dat$Admin_Districts)) /
  length(unique(Admin_Districts$ADM2_EN)) * 100

Percentate_expansion_by_district  ## 42.22222


################################################################################
# 3. SUMMARY STATISTICS — NATIONAL / BUSHENYI / SOROTI
################################################################################

str(TistDist_Dat)

TistDat_summary_stats <- TistDist_Dat %>%
  summarise(
    n_farmers            = n_distinct(Farmer_ID),   # total unique farmers
    n_groups             = n_distinct(Group_ID),     # total unique groups
    n_locations          = n_distinct(Location_ID),  # total unique locations
    n_Districts          = n_distinct(Admin_Districts),
    n_Proj_Area          = n_distinct(Proj_Area),
    n_Village_ID         = n_distinct(Village_ID),
    n_Cluster_ID         = n_distinct(Cluster_ID),
    n_farmers_zero_trees = sum(Trees == 0)           # farmers with 0 trees
  )

TistDat_summary_stats
sum(TistDist_Dat$Trees, na.rm = TRUE)    ## 12573450
sum(TistDist_Dat$Area_Ha, na.rm = TRUE)  ## 26463.25

TistBush_summary_stats <- TistDist_Dat %>%
  filter(Proj_Area == "Bushenyi") %>%
  summarise(
    n_farmers            = n_distinct(Farmer_ID),
    n_groups             = n_distinct(Group_ID),
    n_locations          = n_distinct(Location_ID),
    n_Districts          = n_distinct(Admin_Districts),
    n_Proj_Area          = n_distinct(Proj_Area),
    n_Subcounty          = n_distinct(Subcounty),
    n_Village_ID         = n_distinct(Village_ID),
    n_Cluster_ID         = n_distinct(Cluster_ID),
    n_farmers_zero_trees = sum(Trees == 0)
  )

TistBush_summary_stats
sum(TistDist_Bush$Trees, na.rm = TRUE)    ## 1352838
sum(TistDist_Bush$Area_Ha, na.rm = TRUE)  ## 2947.674

TistSrt_summary_stats <- TistDist_Dat %>%
  filter(Proj_Area == "Soroti") %>%
  summarise(
    n_farmers            = n_distinct(Farmer_ID),
    n_groups             = n_distinct(Group_ID),
    n_locations          = n_distinct(Location_ID),
    n_Districts          = n_distinct(Admin_Districts),
    n_Proj_Area          = n_distinct(Proj_Area),
    n_Subcounty          = n_distinct(Subcounty),
    n_Village_ID         = n_distinct(Village_ID),
    n_Cluster_ID         = n_distinct(Cluster_ID),
    n_farmers_zero_trees = sum(Trees == 0)
  )

TistSrt_summary_stats
sum(TistDist_Soroti$Trees, na.rm = TRUE)    ##
sum(TistDist_Soroti$Area_Ha, na.rm = TRUE)  ##


################################################################################
# 4. SUMMARY STATISTICS TABLE -> MS WORD (.docx)
#    National / Bushenyi / Soroti
#
#    Assumes TistDist_Dat is already loaded in this R session.
#    No file reading happens here.
################################################################################

if (!exists("TistDist_Dat")) {
  stop("Could not find a loaded data object named 'TistDist_Dat' in this session.")
}

## ---- 4.1 Function to calculate summary statistics for one data frame -----

calculate_summary_stats <- function(data) {

  data <- data %>%
    mutate(Years_since_registration = 2024 - reg_date)

  summary_vars <- c(
    "Dist_To_Forest_m",
    "Exposure",
    "Density_winsor99",
    "Trees",
    "Area_Ha",
    "Years_since_registration"
  )

  # preserves display order (chapter order), not alphabetical
  variable_labels <- c(
    Dist_To_Forest_m         = "Distance to forest reserve (m)",
    Exposure                 = "Exposure",
    Density_winsor99         = "Tree density",
    Trees                    = "Tree count",
    Area_Ha                  = "Area (ha)",
    Years_since_registration = "Years since registration"
  )

  missing_vars <- setdiff(summary_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing expected column(s): ", paste(missing_vars, collapse = ", "))
  }

  results <- data %>%
    select(all_of(summary_vars)) %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
    group_by(Variable) %>%
    summarise(
      Mean    = mean(Value, na.rm = TRUE),
      Median  = median(Value, na.rm = TRUE),
      P90     = quantile(Value, probs = 0.90, na.rm = TRUE, names = FALSE),
      Minimum = min(Value, na.rm = TRUE),
      Maximum = max(Value, na.rm = TRUE),
      SD      = sd(Value, na.rm = TRUE),
      N       = sum(!is.na(Value)),
      .groups = "drop"
    ) %>%
    mutate(
      Mean     = round(Mean, 2),
      Median   = round(Median, 2),
      P90      = round(P90, 2),
      SD       = round(SD, 2),
      Range    = paste0(round(Minimum, 2), "\u2013", round(Maximum, 2)),
      Variable = factor(variable_labels[Variable], levels = unname(variable_labels))
    ) %>%
    arrange(Variable) %>%
    select(Variable, Mean, Median, P90, Range, SD, N)

  return(results)
}

## ---- 4.2 Diagnostic: confirm site filters actually match rows ------------

cat("Distinct Admin_Districts values and counts:\n")
print(sort(table(TistDist_Dat$Admin_Districts), decreasing = TRUE))

n_bushenyi <- sum(TistDist_Dat$Admin_Districts == "Bushenyi", na.rm = TRUE)
n_soroti   <- sum(TistDist_Dat$Admin_Districts == "Soroti", na.rm = TRUE)
cat("\nRows matched -- Bushenyi:", n_bushenyi, " | Soroti:", n_soroti, "\n")
cat("(Chapter reports ~2,687 active farmers for Bushenyi and ~2,316 for Soroti,\n",
    "based on the full 12-district / 6-district project areas respectively --\n",
    "compare against those figures before trusting this single-district filter.)\n\n")

## ---- 4.3 Compute summaries -------------------------------------------------

summary_national <- calculate_summary_stats(TistDist_Dat)

summary_bushenyi <- TistDist_Dat %>%
  filter(Admin_Districts == "Bushenyi") %>%
  calculate_summary_stats()

summary_soroti <- TistDist_Dat %>%
  filter(Admin_Districts == "Soroti") %>%
  calculate_summary_stats()

## ---- 4.4 Combine into one long table ---------------------------------------

summary_table <- bind_rows(
  summary_national %>% mutate(Site = "National", .before = 1),
  summary_bushenyi %>% mutate(Site = "Bushenyi", .before = 1),
  summary_soroti   %>% mutate(Site = "Soroti",   .before = 1)
) %>%
  mutate(Site = factor(Site, levels = c("National", "Bushenyi", "Soroti")))

print(summary_table, n = Inf)

## ---- 4.5 Manuscript-style wide table (rows = variables, grouped columns = site) --

summary_table_wide <- summary_table %>%
  select(Site, Variable, Mean, Median, P90, Range, SD) %>%
  pivot_wider(
    names_from  = Site,
    values_from = c(Mean, Median, P90, Range, SD),
    names_glue  = "{Site}_{.value}"
  ) %>%
  select(
    Variable,
    starts_with("National_"),
    starts_with("Bushenyi_"),
    starts_with("Soroti_")
  )

print(summary_table_wide, n = Inf)

## ---- 4.6 Build Word tables (flextable) --------------------------------------
## Three separate tables -- National, Bushenyi, Soroti -- each on its own
## and each fitted to a portrait page.

stats <- c("Mean", "Median", "P90", "Range", "SD")

PAGE_WIDTH <- 6.5  # usable width for a standard portrait page (Letter/A4, 1" margins)

## helper: build a single-site flextable, scaled to fit PAGE_WIDTH
build_site_table <- function(site_df, site_name, caption_text, font_size = 11) {

  flat_names <- c("Variable", paste(site_name, stats, sep = "_"))
  names(site_df) <- flat_names

  ft <- flextable(site_df)

  ft <- set_header_labels(
    ft,
    Variable = "Variable",
    setNames(as.list(stats), flat_names[-1])
  )

  ft <- add_header_row(
    ft,
    top       = TRUE,
    values    = c("", site_name),
    colwidths = c(1, length(stats))
  )

  ft <- merge_h(ft, part = "header")
  ft <- theme_booktabs(ft)
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "header")
  ft <- align(ft, j = 2:ncol(site_df), align = "center", part = "body")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, part = "all")
  ft <- set_caption(ft, caption = caption_text)

  # Give the Variable column a bit more room than the stat columns,
  # then rescale everything proportionally to fit the page width.
  ft <- width(ft, j = 1, width = 1.8)
  ft <- width(ft, j = 2:ncol(site_df), width = (PAGE_WIDTH - 1.8) / (ncol(site_df) - 1))
  ft <- fit_to_width(ft, max_width = PAGE_WIDTH)

  ft
}

## Build each table (Variable + 5 stats = 6 columns each)
ft_national <- build_site_table(
  summary_national %>% select(-N),
  site_name    = "National",
  caption_text = "Table S8a: National summary statistics"
)

ft_bushenyi <- build_site_table(
  summary_bushenyi %>% select(-N),
  site_name    = "Bushenyi",
  caption_text = "Table S8b: Bushenyi summary statistics"
)

ft_soroti <- build_site_table(
  summary_soroti %>% select(-N),
  site_name    = "Soroti",
  caption_text = "Table S8c: Soroti summary statistics"
)

## ---- 4.7 Save all three tables into one Word document (portrait) -----------
## Each table on its own, separated by spacer paragraphs.

doc <- read_docx()
doc <- body_add_flextable(doc, ft_national)
doc <- body_add_par(doc, "")
doc <- body_add_flextable(doc, ft_bushenyi)
doc <- body_add_par(doc, "")
doc <- body_add_flextable(doc, ft_soroti)

print(doc, target = "summary_tables.docx")

cat("\nSaved: summary_tables.docx (Table S8a: National; Table S8b: Bushenyi; Table S8c: Soroti)\n")


################################################################################
# 5. NATIONAL DISTRIBUTION PLOTS — GROUPS, FARMERS, LOCATIONS, VILLAGES
################################################################################

## ---- 5.1 Number of groups per farmer --------------------------------------

TistDat_Farmer_groups <- TistDist_Dat %>%
  group_by(Farmer_ID) %>%
  summarise(n_groups = n_distinct(Group_ID), .groups = "drop")

## Collapse farmers with >=10 groups into one category
TistDat_Farmer_groups_collapsed <- TistDat_Farmer_groups %>%
  mutate(n_groups_cat = ifelse(n_groups >= 10, "\u226510", as.character(n_groups)))

## Summarise number of farmers by number of groups
TistDat_summary_table <- TistDat_Farmer_groups_collapsed %>%
  group_by(n_groups_cat) %>%
  summarise(n_farmers = n(), .groups = "drop") %>%
  mutate(n_groups_cat = factor(n_groups_cat, levels = c(as.character(1:9), "\u226510"))) %>%
  arrange(n_groups_cat)

## Add percentages
TistDat_summary_table_pct <- TistDat_summary_table %>%
  mutate(percent = round(n_farmers / sum(n_farmers) * 100, 1))

TistDat_Farmer_Group_Plot <- ggplot(TistDat_summary_table_pct,
                                     aes(x = n_groups_cat, y = percent)) +
  geom_col(fill = "#1f78b4", width = 0.8) +
  geom_text(aes(label = percent), angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x = "Number of Groups per Farmer",
    y = "Percentage of Farmers"
    # tag = "A"
  ) +
  ylim(0, max(TistDat_summary_table_pct$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

TistDat_Farmer_Group_Plot

## ---- 5.2 Number of farmers per group ---------------------------------------

TistDat_Group_farmers <- TistDist_Dat %>%
  group_by(Group_ID) %>%
  summarise(n_farmers = n_distinct(Farmer_ID), .groups = "drop")

## Collapse groups: 1-12 individually, >12 as one bin
TistDat_Group_farmers_collapsed <- TistDat_Group_farmers %>%
  mutate(n_farmers_cat = case_when(
    n_farmers <= 12 ~ as.character(n_farmers),
    n_farmers > 12  ~ "\u226513"
  ))

## Summarise number of groups by (collapsed) farmer count
TistDat_summary_group_farmers <- TistDat_Group_farmers_collapsed %>%
  group_by(n_farmers_cat) %>%
  summarise(n_groups = n(), .groups = "drop") %>%
  mutate(n_farmers_cat = factor(n_farmers_cat, levels = c(as.character(1:12), "\u226513"))) %>%
  arrange(n_farmers_cat)

## Add percentages
TistDat_summary_group_farmers_pct <- TistDat_summary_group_farmers %>%
  mutate(percent = round(n_groups / sum(n_groups) * 100, 1))
TistDat_summary_group_farmers_pct

TistDat_Group_Farmer_Plot <- ggplot(TistDat_summary_group_farmers_pct,
                                     aes(x = n_farmers_cat, y = percent)) +
  geom_col(fill = "#1f78b4", width = 0.8) +
  geom_text(aes(label = percent), angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x = "Number of Farmers per Group",
    y = "Percentage of Groups"
  ) +
  ylim(0, max(TistDat_summary_group_farmers_pct$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

TistDat_Group_Farmer_Plot

## ---- 5.3 Number of grove locations per farmer -------------------------------

TistDat_Farmer_locations <- TistDist_Dat %>%
  group_by(Farmer_ID) %>%
  summarise(n_locations = n_distinct(Location_ID), .groups = "drop")

TistDat_locations_summary_binned <- TistDat_Farmer_locations %>%
  mutate(location_bin = ifelse(n_locations >= 10, "\u226510", as.character(n_locations))) %>%
  group_by(location_bin) %>%
  summarise(n_farmers = n(), .groups = "drop") %>%
  mutate(
    percent      = round(n_farmers / sum(n_farmers) * 100, 1),
    location_bin = factor(location_bin, levels = c(as.character(1:9), "\u226510"))
  )

TistDat_Farmer_Location_Plot <- ggplot(TistDat_locations_summary_binned,
                                        aes(x = location_bin, y = percent)) +
  geom_col(fill = "darkgreen", width = 0.8) +
  geom_text(aes(label = percent), angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x = "Number of Grove Locations per Farmer",
    y = "Percentage of Farmers"
    # tag = "C"
  ) +
  ylim(0, max(TistDat_locations_summary_binned$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

TistDat_Farmer_Location_Plot

## ---- 5.4 Number of farmer groves per village ---------------------------------

TistDat_Location_farmers <- TistDist_Dat %>%
  group_by(Village_ID) %>%
  summarise(n_farmers = n_distinct(Location_ID), .groups = "drop")

TistDat_summary_location_farmers <- TistDat_Location_farmers %>%
  group_by(n_farmers) %>%
  summarise(n_locations = n(), .groups = "drop") %>%
  arrange(n_farmers)

TistDat_summary_location_farmers_binned <- TistDat_Location_farmers %>%
  mutate(
    farmer_bin = case_when(
      n_farmers %in% 1:29 ~ as.character(n_farmers),
      n_farmers > 29      ~ "\u226530"
    )
  ) %>%
  group_by(farmer_bin) %>%
  summarise(n_locations = n(), .groups = "drop") %>%
  mutate(
    percent    = round(n_locations / sum(n_locations) * 100, 1),
    farmer_bin = factor(farmer_bin, levels = c(as.character(1:30), "\u226530"))
  )

TistDat_Location_Farmer_Plot <- ggplot(TistDat_summary_location_farmers_binned,
                                        aes(x = farmer_bin, y = percent)) +
  geom_col(fill = "darkgreen", width = 0.8) +
  geom_text(aes(label = percent), angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x = "Number of Farmer Groves per Village",
    y = "Percentage of Villages"
    # tag = "D"
  ) +
  ylim(0, max(TistDat_summary_location_farmers_binned$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

TistDat_Location_Farmer_Plot

## ---- 5.5 Combined 4-panel national figure ------------------------------------

TistDat_Combined_Plot <- wrap_plots(
  TistDat_Farmer_Group_Plot, TistDat_Group_Farmer_Plot,
  TistDat_Farmer_Location_Plot, TistDat_Location_Farmer_Plot,
  ncol = 2, nrow = 2,
  widths = c(1, 2)  # push the ratio further, since 1.6 wasn't visually distinct
) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag          = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
  )

TistDat_Combined_Plot

ggsave(
  filename = "Output/Manuscript 3 graphs/TistDat_Combined_Plot.png",
  plot     = TistDat_Combined_Plot,
  width    = 200,  # mm — double-column width, same as before
  height   = 200,  # mm — roughly doubled from the 2-panel version to accommodate 2 rows
  units    = "mm",
  dpi      = 600
)


################################################################################
# 6. BUSHENYI vs SOROTI DISTRIBUTION PLOTS — GROUPS, FARMERS, LOCATIONS
################################################################################

## ---- 6.1 Number of groups per farmer, by project area ------------------------

Farmer_groups <- Tist_Bush_Soroti %>%
  group_by(Farmer_ID, Proj_Area) %>%
  summarise(n_groups = n_distinct(Group_ID), .groups = "drop")

Farmer_groups_collapsed <- Farmer_groups %>%
  mutate(n_groups_cat = ifelse(n_groups >= 10, "\u226510", as.character(n_groups)))

summary_table <- Farmer_groups_collapsed %>%
  group_by(Proj_Area, n_groups_cat) %>%
  summarise(n_farmers = n(), .groups = "drop") %>%
  mutate(n_groups_cat = factor(n_groups_cat, levels = c(as.character(1:9), "\u226510"))) %>%
  arrange(Proj_Area, n_groups_cat)
summary_table

summary_table_pct <- summary_table %>%
  group_by(Proj_Area) %>%
  mutate(percent = round(n_farmers / sum(n_farmers) * 100, 1))
summary_table_pct

Farmer_Group_Plot <- ggplot(summary_table_pct,
                             aes(x = n_groups_cat, y = percent, fill = Proj_Area)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = percent), position = position_dodge(width = 0.9),
            angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x    = "Number of Groups per Farmer",
    y    = "Percentage of Farmers",
    fill = "Project Area"
    # tag = "A"
  ) +
  scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
  ylim(0, max(summary_table_pct$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    legend.position   = "bottom",
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 12),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Farmer_Group_Plot

## ---- 6.2 Number of farmers per group, by project area ------------------------

Group_farmers <- Tist_Bush_Soroti %>%
  group_by(Group_ID, Proj_Area) %>%
  summarise(n_farmers = n_distinct(Farmer_ID), .groups = "drop")

## Collapse groups: 1-12 individually, >12 as one bin
Group_farmers_collapsed <- Group_farmers %>%
  mutate(n_farmers_cat = case_when(
    n_farmers <= 12 ~ as.character(n_farmers),
    n_farmers > 12  ~ "\u226513"
  ))

summary_group_farmers <- Group_farmers_collapsed %>%
  group_by(Proj_Area, n_farmers_cat) %>%
  summarise(n_groups = n(), .groups = "drop") %>%
  mutate(n_farmers_cat = factor(n_farmers_cat, levels = c(as.character(1:12), "\u226513"))) %>%
  arrange(Proj_Area, n_farmers_cat)
summary_group_farmers

summary_group_farmers_pct <- summary_group_farmers %>%
  group_by(Proj_Area) %>%
  mutate(percent = round(n_groups / sum(n_groups) * 100, 1))
summary_group_farmers_pct

Group_Farmer_Pct_Plot <- ggplot(summary_group_farmers_pct,
                                 aes(x = n_farmers_cat, y = percent, fill = Proj_Area)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = percent), position = position_dodge(width = 0.9),
            angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x    = "Number of Farmers per Group",
    y    = "Percentage of Groups",
    fill = "Project Area"
    # tag = "B"
  ) +
  scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
  ylim(0, max(summary_group_farmers_pct$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    legend.position   = "bottom",
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 12),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Group_Farmer_Pct_Plot

## ---- 6.3 Number of groves per village, by project area -----------------------

Location_farmers <- Tist_Bush_Soroti %>%
  group_by(Village_ID, Proj_Area) %>%
  summarise(n_farmers = n_distinct(Location_ID), .groups = "drop")

summary_location_farmers <- Location_farmers %>%
  group_by(Proj_Area, n_farmers) %>%
  summarise(n_locations = n(), .groups = "drop") %>%
  arrange(Proj_Area, n_farmers)

summary_location_farmers_pct <- summary_location_farmers %>%
  group_by(Proj_Area) %>%
  mutate(percent = round(n_locations / sum(n_locations) * 100, 1))

summary_location_farmers_binned <- summary_location_farmers_pct %>%
  mutate(
    n_farmers_bin = case_when(
      n_farmers < 20  ~ as.character(n_farmers),
      n_farmers >= 20 ~ "\u226520"
    )
  )

## Re-aggregate percentages by binned values
summary_location_farmers_binned_pct <- summary_location_farmers_binned %>%
  group_by(Proj_Area, n_farmers_bin) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  mutate(n_farmers_bin = factor(n_farmers_bin, levels = c(as.character(1:19), "\u226520")))

Location_Farmer_Pct_Plot <- ggplot(summary_location_farmers_binned_pct,
                                    aes(x = n_farmers_bin, y = percent, fill = Proj_Area)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = percent), position = position_dodge(width = 0.9),
            angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x    = "Number of Groves per Village",
    y    = "Percentage of Villages",
    fill = "Project Area"
  ) +
  scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
  ylim(0, max(summary_location_farmers_binned_pct$percent) * 1.3) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    legend.position   = "bottom",
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 12),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Location_Farmer_Pct_Plot

## ---- 6.4 Number of groves per farmer, by project area ------------------------

Farmer_locations <- Tist_Bush_Soroti %>%
  group_by(Farmer_ID, Proj_Area) %>%
  summarise(n_locations = n_distinct(Location_ID), .groups = "drop")

## Collapse farmers with >=10 locations into one category, consistent with
## single-region Panel C
Farmer_locations_collapsed <- Farmer_locations %>%
  mutate(n_locations_cat = ifelse(n_locations >= 10, "\u226510", as.character(n_locations)))

summary_farmer_locations <- Farmer_locations_collapsed %>%
  group_by(Proj_Area, n_locations_cat) %>%
  summarise(n_farmers = n(), .groups = "drop") %>%
  mutate(n_locations_cat = factor(n_locations_cat, levels = c(as.character(1:9), "\u226510"))) %>%
  arrange(Proj_Area, n_locations_cat)

summary_farmer_locations_pct <- summary_farmer_locations %>%
  group_by(Proj_Area) %>%
  mutate(percent = round(n_farmers / sum(n_farmers) * 100, 1))

Farmer_Location_Pct_Plot <- ggplot(summary_farmer_locations_pct,
                                    aes(x = n_locations_cat, y = percent, fill = Proj_Area)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = percent), position = position_dodge(width = 0.9),
            angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x    = "Number of Groves per Farmer",
    y    = "Percentage of Farmers",
    fill = "Project Area"
    # tag = "D"
  ) +
  scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
  ylim(0, max(summary_farmer_locations_pct$percent) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 12),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    legend.position   = "bottom",
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 12),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Farmer_Location_Pct_Plot

## ---- 6.5 Combined 4-panel Bushenyi/Soroti figure ------------------------------
## All four plots already have their tag = "A"/"B"/"C"/"D" commented out in
## the document provided, so plot_annotation() below is the single source
## of truth for panel lettering — no further edits needed to labs() calls.

BushSrt_Combined_Plot <- wrap_plots(
  Farmer_Group_Plot, Group_Farmer_Pct_Plot,
  Farmer_Location_Pct_Plot, Location_Farmer_Pct_Plot,
  ncol = 2, nrow = 2,
  widths = c(1, 2)  # push the ratio further, since 1.6 wasn't visually distinct
) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag          = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
  )

BushSrt_Combined_Plot

ggsave(
  filename = "Output/Manuscript 3 graphs/BushSrt_Combined_Plot.png",
  plot     = BushSrt_Combined_Plot,
  width    = 200,  # mm — double-column width, matches other combined figures
  height   = 200,  # mm — doubled from the 2-panel version to accommodate 2 rows
  units    = "mm",
  dpi      = 600
)


################################################################################
# 7. SUBCOUNTY PENETRATION & CLUSTER COMPOSITION
################################################################################

## ---- 7.1 Proportion of population who joined the program, by subcounty ------

BushSoroti_TISTsubcountypopn <- read.csv("Data/TISTDat/BushSoroti_cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered_communicationoptions.csv") %>%
  distinct()

farmer_per_subcounty <- BushSoroti_TISTsubcountypopn %>%
  group_by(Proj_Area, Subcounty) %>%
  summarise(
    n_farmers  = n_distinct(Farmer_ID),
    population = first(Total_Subcounty_popn.),
    .groups = "drop"
  ) %>%
  mutate(farmer_percent = (n_farmers / population) * 100) %>%
  arrange(desc(farmer_percent))
farmer_per_subcounty

summary(farmer_per_subcounty$farmer_percent)

## ---- Helper: penetration bar plot, with small-category collapsing -----------

plot_penetration <- function(data, area_name,
                              threshold  = NULL,
                              fill_color = "#1f78b4",
                              tag        = NULL) {

  plot_data <- data %>%
    filter(Proj_Area == area_name)

  ## Collapse small categories if threshold supplied
  if (!is.null(threshold)) {

    plot_data <- plot_data %>%
      mutate(
        Subcounty_grouped = ifelse(
          farmer_percent < threshold,
          paste0("Other (<", threshold, "%)"),
          Subcounty
        )
      ) %>%
      group_by(Subcounty_grouped) %>%
      summarise(
        n_farmers  = sum(n_farmers),
        population = sum(population),
        .groups = "drop"
      ) %>%
      mutate(farmer_percent = (n_farmers / population) * 100)

  } else {

    plot_data <- plot_data %>%
      mutate(Subcounty_grouped = Subcounty)
  }

  ## Order bars, forcing "Other" to be last regardless of its value
  other_label <- paste0("Other (<", threshold, "%)")

  plot_data <- plot_data %>%
    arrange(desc(farmer_percent)) %>%
    mutate(
      Subcounty_grouped = factor(
        Subcounty_grouped,
        levels = c(setdiff(Subcounty_grouped, other_label), other_label)
      )
    )

  ggplot(plot_data, aes(x = Subcounty_grouped, y = farmer_percent)) +
    geom_col(fill = fill_color, width = 0.7) +
    geom_text(aes(label = round(farmer_percent, 2)),
              angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
    labs(
      x   = "Sub-county",
      y   = "% of Population in TIST",
      tag = tag
    ) +
    expand_limits(y = max(plot_data$farmer_percent) * 1.15) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin       = margin(t = 15, l = 5, r = 5, b = 5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.border      = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.text         = element_text(size = 10),
      axis.text.x       = element_text(angle = 90, hjust = 1),
      axis.title        = element_text(size = 13),
      plot.tag          = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.02, 0.98)
    )
}

Sorotipenetration_plot <- plot_penetration(
  farmer_per_subcounty,
  area_name  = "Soroti",
  threshold  = 0.1,
  fill_color = "#33a02c"
  # tag        = "A"
)
Sorotipenetration_plot

Bushenyipenetration_plot <- plot_penetration(
  farmer_per_subcounty,
  area_name  = "Bushenyi",
  threshold  = 0.5,
  fill_color = "#1f78b4"
  # tag        = "B"
)
Bushenyipenetration_plot

## ---- 7.2 Cluster composition -------------------------------------------------

cluster_summary <- Tist_Bush_Soroti %>%
  group_by(Proj_Area, Cluster_ID) %>%
  summarise(n_farmers = n_distinct(Farmer_ID), .groups = "drop") %>%
  group_by(Proj_Area) %>%
  mutate(percent_of_farmers = (n_farmers / sum(n_farmers)) * 100) %>%
  ungroup() %>%
  arrange(Proj_Area, desc(n_farmers))

head(cluster_summary)
summary(cluster_summary$n_farmers)
summary(cluster_summary$percent_of_farmers)

## ---- Soroti clusters ----------------------------------------------------------
Soroti_cluster_plot_data <- cluster_summary %>%
  filter(Proj_Area == "Soroti") %>%
  mutate(Cluster_grouped = ifelse(percent_of_farmers < 1, "Other (<1%)", Cluster_ID)) %>%
  group_by(Cluster_grouped) %>%
  summarise(
    n_farmers          = sum(n_farmers),
    percent_of_farmers = sum(percent_of_farmers),
    .groups = "drop"
  ) %>%
  arrange(desc(percent_of_farmers)) %>%
  mutate(Cluster_grouped = factor(Cluster_grouped, levels = Cluster_grouped))

Soroti_cluster_plot <- ggplot(Soroti_cluster_plot_data,
                               aes(x = Cluster_grouped, y = percent_of_farmers)) +
  geom_col(fill = "#33a02c", width = 0.7) +
  geom_text(aes(label = round(percent_of_farmers, 2)),
            angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x = "Cluster",
    y = "% of TIST Farmers"
    # tag = "C"   # set to whatever letter fits your panel sequence
  ) +
  expand_limits(y = max(Soroti_cluster_plot_data$percent_of_farmers) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 10),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Soroti_cluster_plot

## ---- Bushenyi clusters ---------------------------------------------------------
Bushenyi_cluster_plot_data <- cluster_summary %>%
  filter(Proj_Area == "Bushenyi") %>%
  mutate(Cluster_grouped = ifelse(percent_of_farmers < 2, "Other (<2%)", Cluster_ID)) %>%
  group_by(Cluster_grouped) %>%
  summarise(
    n_farmers          = sum(n_farmers),
    percent_of_farmers = sum(percent_of_farmers),
    .groups = "drop"
  ) %>%
  arrange(desc(percent_of_farmers)) %>%
  ## Force "Other" to be last on x-axis
  mutate(
    Cluster_grouped = factor(
      Cluster_grouped,
      levels = c(setdiff(Cluster_grouped, "Other (<2%)"), "Other (<2%)")
    )
  )

Bushenyi_cluster_plot <- ggplot(Bushenyi_cluster_plot_data,
                                 aes(x = Cluster_grouped, y = percent_of_farmers)) +
  geom_col(fill = "#1f78b4", width = 0.7) +
  geom_text(aes(label = round(percent_of_farmers, 2)),
            angle = 90, hjust = -0.1, vjust = 0.5, size = 3) +
  labs(
    x = "Cluster",
    y = "% of TIST Farmers"
    # tag = "D"
  ) +
  expand_limits(y = max(Bushenyi_cluster_plot_data$percent_of_farmers) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 10),
    axis.text.x       = element_text(angle = 90, hjust = 1),
    axis.title        = element_text(size = 13),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Bushenyi_cluster_plot

## ---- 7.3 Combined 4-panel penetration + cluster figure ------------------------
## Reminder: tag = ... is already commented out in all four individual plot
## objects (Sorotipenetration_plot, Bushenyipenetration_plot, Soroti_cluster_plot,
## Bushenyi_cluster_plot), so plot_annotation() below is the single source of
## truth for panel lettering — no further edits needed.

Penetration_Cluster_Combined_Plot <-
  Sorotipenetration_plot + Bushenyipenetration_plot +
  Soroti_cluster_plot + Bushenyi_cluster_plot +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "A") &     # capital letters, matches convention used elsewhere
  theme(
    plot.margin       = margin(t = 15, l = 5, r = 5, b = 5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 1.05)     # upper-left corner, consistent with A-D panel convention
  )

Penetration_Cluster_Combined_Plot

ggsave(
  filename = "Output/Manuscript 3 graphs/Penetration_Cluster_Combined_Plot.png",
  plot     = Penetration_Cluster_Combined_Plot,
  width    = 200,  # mm — double-column width, matches other combined figures
  height   = 200,  # mm — matches other 4-panel figures
  units    = "mm",
  dpi      = 600
)


################################################################################
# 8. GROUPS BY NUMBER OF LOCATIONS (BUSHENYI vs SOROTI)
################################################################################

Group_locations <- Tist_Bush_Soroti %>%
  group_by(Group_ID, Proj_Area) %>%
  summarise(n_locations = n_distinct(Location_ID), .groups = "drop")

summary_group_locations <- Group_locations %>%
  group_by(Proj_Area, n_locations) %>%
  summarise(n_groups = n(), .groups = "drop") %>%
  arrange(Proj_Area, n_locations)

summary_group_locations

summary_group_locations_pct <- summary_group_locations %>%
  group_by(Proj_Area) %>%
  mutate(percent = round(n_groups / sum(n_groups) * 100, 1))

Group_Location_pct_plot <- ggplot(summary_group_locations_pct,
                                   aes(x = factor(n_locations), y = percent, fill = Proj_Area)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = percent), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 4) +
  labs(
    x     = "Number of Locations per Group",
    y     = "Percentage of Groups",
    fill  = "Project Area",
    title = "Percentage of Groups by Locations"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
  ylim(0, max(summary_group_locations_pct$percent) * 1.15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 13),
    legend.position  = "bottom"
  )

Group_Location_pct_plot

ggsave("Output/Manuscript 3 graphs/BushSrtGroup_Location_pct_plot.png",
       Group_Location_pct_plot, dpi = 600, width = 12, height = 6)


################################################################################
# 9. FIGURE S1 (SUPPLEMENTARY) — COMMUNICATION CHANNELS: BUSHENYI vs SOROTI
################################################################################
#
# Reload the communications-specific dataset separately here, rather than
# carrying its extra ~70 neighbour-threshold columns through the rest of
# the script. Only Farmer_ID, Proj_Area, and the communication-channel
# columns are needed for this figure.
################################################################################

TistDist_Comm <- read.csv("Data/TISTDat/BushSoroti_cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered_communicationoptions.csv") %>%
  rename(Admin_Districts = Admin_Districts.x) %>%
  select(-Admin_Districts.y) %>%
  distinct()

## Aggregate to one row per Subcounty (communication figures are already
## subcounty-level census counts, duplicated across every farmer row in the
## merge) before computing site-level shares, to avoid double-counting.

comm_vars <- c(
  "Radio", "Word_of_Mouth", "Phone_calls", "TV",
  "Community_meetings", "Social_media",
  "Community_Announcer", "Print_Media",
  "Other_sources"
)

Comm_subcounty <- TistDist_Comm %>%
  filter(Proj_Area %in% c("Bushenyi", "Soroti")) %>%
  distinct(Proj_Area, Subcounty, Total_Households, across(all_of(comm_vars))) %>%
  mutate(across(all_of(comm_vars), ~ as.numeric(gsub(",", "", as.character(.)))))

head(Comm_subcounty)

Comm_long <- Comm_subcounty %>%
  pivot_longer(cols = all_of(comm_vars), names_to = "Channel", values_to = "n_households")

Comm_summary <- Comm_long %>%
  group_by(Proj_Area, Channel) %>%
  summarise(
    total_channel_households = sum(n_households, na.rm = TRUE),
    total_households          = sum(Total_Households, na.rm = TRUE) / n_distinct(Channel),
    mean_pct                  = 100 * total_channel_households / total_households,
    .groups = "drop"
  )

channel_order <- Comm_summary %>%
  group_by(Channel) %>%
  summarise(overall = mean(mean_pct), .groups = "drop") %>%
  arrange(desc(overall)) %>%
  pull(Channel)

Comm_summary <- Comm_summary %>%
  mutate(Channel = factor(Channel, levels = channel_order))

FigureS1_Communication_Plot <- ggplot(Comm_summary,
                                       aes(x = Channel, y = mean_pct, fill = Proj_Area)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = round(mean_pct, 1)), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3) +
  labs(
    x    = "Communication channel",
    y    = "% of households",
    fill = "Project Area"
  ) +
  scale_fill_manual(values = c("Bushenyi" = "#1f78b4", "Soroti" = "#33a02c")) +
  ylim(0, max(Comm_summary$mean_pct) * 1.15) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title       = element_text(size = 13),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 12)
  )

FigureS1_Communication_Plot

ggsave(
  filename = "C:/workspace/Emenyu_et_al_TIST_leveraging_adoption/Output/Manuscript 3 graphs/manuscript_new_plots/FigureS1_Communication_Plot.png",
  plot     = FigureS1_Communication_Plot,
  width    = 183,
  height   = 120,
  units    = "mm",
  dpi      = 600
)


################################################################################
# 10. SECTION 4.2 — RIDGELINE PLOTS: TREE COUNT, AREA, DENSITY
################################################################################
#
# Uses TistDist_Dat (already loaded upstream) directly — no upstream
# load-block changes needed.
################################################################################

TistDist_Dat_plot <- TistDist_Dat %>%
  mutate(
    Site = case_when(
      Proj_Area %in% c("Bushenyi", "Soroti") ~ Proj_Area,
      TRUE ~ "National"
    )
  ) %>%
  mutate(Site = factor(Site, levels = c("National", "Bushenyi", "Soroti")))

## ---- 10.1 Panel A: Tree count -------------------------------------------------

Trees_plot_dat <- TistDist_Dat_plot %>%
  filter(!is.na(Trees), Trees > 0)

stats_trees <- Trees_plot_dat %>%
  group_by(Site) %>%
  summarise(
    median_val = median(Trees),
    mean_val   = mean(Trees),
    p90_val    = quantile(Trees, 0.9),
    .groups = "drop"
  )

Ridgeline_Tree_Plot <- ggplot(Trees_plot_dat, aes(x = Trees, y = Site, fill = Site)) +
  geom_density_ridges(
    scale = 1.2, alpha = 0.7, quantile_lines = TRUE,
    quantiles = c(0.5, 0.9), rel_min_height = 0.01
  ) +
  geom_point(data = stats_trees, aes(x = mean_val, y = Site),
             inherit.aes = FALSE, size = 2, color = "black") +
  geom_text(data = stats_trees,
            aes(x = median_val, y = Site, label = paste0("Med: ", median_val)),
            inherit.aes = FALSE, vjust = -0.5, hjust = 0, size = 3) +
  geom_text(data = stats_trees,
            aes(x = mean_val, y = Site, label = paste0("Mean: ", round(mean_val, 1))),
            inherit.aes = FALSE, vjust = 1.5, hjust = 0, size = 3) +
  geom_text(data = stats_trees,
            aes(x = p90_val, y = Site, label = paste0("P90: ", round(p90_val, 0))),
            inherit.aes = FALSE, vjust = -1.5, hjust = 0, size = 3) +
  scale_x_log10() +
  theme_ridges(font_size = 11) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  labs(x = "Number of trees per farmer (log scale)")

## ---- 10.2 Panel B: Area --------------------------------------------------------

Area_plot_dat <- TistDist_Dat_plot %>%
  filter(!is.na(Area_Ha), Area_Ha > 0)

min_area <- min(Area_plot_dat$Area_Ha)

stats_area <- Area_plot_dat %>%
  group_by(Site) %>%
  summarise(
    median_val = median(Area_Ha),
    mean_val   = mean(Area_Ha),
    p90_val    = quantile(Area_Ha, 0.9),
    .groups = "drop"
  )

Ridgeline_Area_Plot <- ggplot(Area_plot_dat, aes(x = Area_Ha, y = Site, fill = Site)) +
  geom_density_ridges(
    scale = 1.2, alpha = 0.7, quantile_lines = TRUE,
    quantiles = c(0.5, 0.9), rel_min_height = 0.01
  ) +
  geom_point(data = stats_area, aes(x = mean_val, y = Site),
             inherit.aes = FALSE, size = 2, color = "black") +
  geom_text(data = stats_area,
            aes(x = median_val, y = Site, label = paste0("Med: ", round(median_val, 2))),
            inherit.aes = FALSE, vjust = -0.5, hjust = 0, size = 3) +
  geom_text(data = stats_area,
            aes(x = mean_val, y = Site, label = paste0("Mean: ", round(mean_val, 2))),
            inherit.aes = FALSE, vjust = 1.5, hjust = 0, size = 3) +
  geom_text(data = stats_area,
            aes(x = p90_val, y = Site, label = paste0("P90: ", round(p90_val, 2))),
            inherit.aes = FALSE, vjust = -1.5, hjust = 0, size = 3) +
  scale_x_log10(limits = c(min_area, NA)) +
  theme_ridges(font_size = 11) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  labs(x = "Farm area (ha, log scale)")

## ---- 10.3 Panel C: Density ------------------------------------------------------

Density_plot_dat <- TistDist_Dat_plot %>%
  filter(!is.na(Density_winsor99), Density_winsor99 > 0)

stats_density <- Density_plot_dat %>%
  group_by(Site) %>%
  summarise(
    median_val = median(Density_winsor99),
    mean_val   = mean(Density_winsor99),
    p90_val    = quantile(Density_winsor99, 0.9),
    .groups = "drop"
  )

min_density <- min(Density_plot_dat$Density_winsor99, na.rm = TRUE)

Ridgeline_Density_Plot <- ggplot(Density_plot_dat, aes(x = Density_winsor99, y = Site, fill = Site)) +
  geom_density_ridges(
    scale = 1.2, alpha = 0.7, quantile_lines = TRUE,
    quantiles = c(0.5, 0.9), rel_min_height = 0.01,
    trim = TRUE, from = min_density
  ) +
  geom_point(data = stats_density, aes(x = mean_val, y = Site),
             inherit.aes = FALSE, size = 2, color = "black") +
  geom_text(data = stats_density,
            aes(x = median_val, y = Site, label = paste0("Med: ", round(median_val, 0))),
            inherit.aes = FALSE, vjust = -0.5, hjust = 0, size = 3) +
  geom_text(data = stats_density,
            aes(x = mean_val, y = Site, label = paste0("Mean: ", round(mean_val, 0))),
            inherit.aes = FALSE, vjust = 1.5, hjust = 0, size = 3) +
  geom_text(data = stats_density,
            aes(x = p90_val, y = Site, label = paste0("P90: ", round(p90_val, 0))),
            inherit.aes = FALSE, vjust = -1.5, hjust = 0, size = 3) +
  scale_x_log10() +
  theme_ridges(font_size = 11) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  labs(x = "Planted tree density, 99th percentile capped (trees per hectare, log scale)")

## ---- 10.4 Combine into single 3-panel figure, tagged A, B, C -------------------

Figure5_Ridgeline_Plot <-
  Ridgeline_Tree_Plot + Ridgeline_Area_Plot + Ridgeline_Density_Plot +
  plot_layout(ncol = 1, nrow = 3) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

Figure5_Ridgeline_Plot

ggsave(
  filename = "C:/workspace/Emenyu_et_al_TIST_leveraging_adoption/Output/Manuscript 3 graphs/manuscript_new_plots/Figure5_Ridgeline_Plot.png",
  plot     = Figure5_Ridgeline_Plot,
  width    = 200,
  height   = 260,
  units    = "mm",
  dpi      = 600
)

################################################################################
# END OF SCRIPT
################################################################################
