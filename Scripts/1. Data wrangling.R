################################################################################
# TIST Uganda — Data Cleaning, pseudonymisation, and Spatial Feature Construction
################################################################################

## ---- Packages ---------------------------------------------------------- ##
library(Rglpk)
library(ergm)
library(tidyverse)
library(sf)
library(geosphere)
library(stringr)
library(statnet)
library(tergm)
library(ergm.multi)
library(network)
library(ggcorrplot)
library(ggpubr)
library(patchwork)
library(modelsummary)
library(dbscan)
library(tibble)


################################################################################
# 1. LOAD RAW DATA
################################################################################

TistDat4 <- read.csv("C:/Users/ae474/OneDrive - University of Exeter/PhD work in progress/Data/R Large Data/tistDat4.csv") %>%
  drop_na() %>%
  distinct()

str(TistDat4)

## Check number of districts per region
TistDat_Reg_DistCounts <- TistDat4 %>%
  select(Region, Admin.Districts) %>%
  na.omit() %>%
  distinct() %>%
  arrange(Region) %>%
  group_by(Region) %>%
  summarise(n_districts = n(), .groups = "drop")


################################################################################
# 2. NAME CLEANING & DUPLICATE IDENTIFICATION
################################################################################

TistDat5 <- TistDat4 %>%
  mutate(
    name_original = Name,
    name_clean = Name %>%
      str_to_lower() %>%
      str_replace_all("[\\.,]", "") %>%
      str_squish() %>%
      str_remove("\\b(\\d+|[a-z]{1}|[ivxlc]{1,4})\\s*$") %>%  # remove trailing suffix
      str_remove("\\d+$") %>%                                  # remove digits directly attached (e.g. "stephen1")
      str_squish()
  )

str(TistDat5)
colnames(TistDat5)

## Potential duplicate records
potential_duplicates <- TistDat5 %>%
  select("Name", "name_clean", "Visual.Description", "Grove.Area.in.Ha", "Trees",
         "TIST.Number", "Group.Name", "Cluster", "Village", "Group.Center", "Area",
         "reg_date", "Admin.Districts", "Region", "Latitude_dd", "Longitude_dd") %>%
  group_by(name_clean, Cluster, Group.Center, Admin.Districts) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(name_clean)

## Select and rename core columns
TistDat6 <- TistDat5 %>%
  select("name_clean", "Visual.Description", "Grove.Area.in.Ha", "Trees", "TIST.Number",
         "Group.Name", "Cluster", "Village", "Group.Center", "Area", "reg_date",
         "Admin.Districts", "Region", "Latitude_dd", "Longitude_dd") %>%
  set_names(c("Name", "Visual_Description", "Area_Ha", "Trees", "TIST_No", "Group_No",
              "Cluster", "Village", "Group_Center", "Proj_Area", "reg_date",
              "Admin_Districts", "Region", "latitude", "longitude")) %>%
  distinct()

TistDat6 <- TistDat6 %>%
  rename(
    Group_Name = Group_No,   # this had names
    Group_No   = TIST_No     # this is the actual group number
  )

TistDat6 %>%
  summarise(
    empty_names   = sum(Name == "", na.rm = TRUE),
    missing_names = sum(is.na(Name))
  )  ## 4 names


################################################################################
# 3. MERGE MULTI-GROVE FARMERS (GROVES WITHIN 50M TREATED AS ONE GROVE)
################################################################################

## Some farmers have multiple groves at the same site (e.g. by the roadside,
## by the house). To avoid misrepresentation, groves within 50m of each other
## are merged into a single grove.

## Identify farmers with more than one grove
multi_grove_farmers <- TistDat6 %>%
  group_by(Name) %>%
  filter(n() > 1) %>%
  ungroup()

## Check how many groves each farmer has
multi_grove_farmers %>%
  count(Name) %>%
  arrange(desc(n))

## Identify whose multiple groves fall within 50m, 75m, 100m, 150m, 200m
mg_sf <- st_as_sf(
  multi_grove_farmers,
  coords = c("longitude", "latitude"),
  crs = 4326
) %>%
  st_transform(32736)

## Helper: flag whether a farmer has any two groves within a given distance
check_within_farmer <- function(sf_obj, dist_m) {
  dist_u <- units::set_units(dist_m, "m")
  zero_u <- units::set_units(0, "m")

  sf_obj %>%
    group_by(Name) %>%
    mutate(
      !!paste0("within_", dist_m, "m") := {
        d <- st_distance(geometry)
        any(d > zero_u & d <= dist_u)
      }
    ) %>%
    ungroup()
}

mg_sf <- mg_sf %>%
  check_within_farmer(50) %>%
  check_within_farmer(75) %>%
  check_within_farmer(100) %>%
  check_within_farmer(150) %>%
  check_within_farmer(200)

## ---- Summary and plot: share of farmers with clustered groves ---------- ##
farmer_summary <- mg_sf %>%
  st_drop_geometry() %>%
  distinct(Name, within_50m, within_75m, within_100m, within_150m, within_200m)

N_total <- nrow(farmer_summary)  # total number of multi-grove farmers

dist_summary <- tibble(
  Distance = c("50 m", "75 m", "100 m", "150 m", "200 m"),
  Count = c(
    sum(farmer_summary$within_50m),
    sum(farmer_summary$within_75m),
    sum(farmer_summary$within_100m),
    sum(farmer_summary$within_150m),
    sum(farmer_summary$within_200m)
  )
) %>%
  mutate(Percent = 100 * Count / N_total) %>%
  mutate(
    Distance = factor(
      Distance,
      levels = c("50 m", "75 m", "100 m", "150 m", "200 m", "> 200 m")
    )
  )

dist_summary

farmergrove_distances <- ggplot(dist_summary, aes(x = Distance, y = Percent)) +
  geom_col() +
  geom_text(
    aes(label = paste0(round(Percent, 1), "% (n=", Count, ")")),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Multi-grove farmers with clustered groves",
    subtitle = paste0("Percent of farmers with \u22652 groves within distance (N = ", N_total, ")"),
    y = "Percent of farmers",
    x = "Distance threshold"
  ) +
  ylim(0, max(dist_summary$Percent) * 1.15) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(size = 10),
    axis.text.x      = element_text(angle = 90, hjust = 1),
    axis.title       = element_text(size = 13)
  )

farmergrove_distances  ## 50m (10.6%), 75m (19.4%), 100m (29%), 150m (41.4%), 200m (47.9%)

ggsave("Output/Manuscript 3 graphs/farmergrove_distances.png",
       farmergrove_distances, dpi = 600, width = 12, height = 6)

## ---- Assign placeholder names for missing entries ----------------------- ##
TistDat6 <- TistDat6 %>%
  mutate(
    Name = ifelse(Name == "" | is.na(Name),
                  paste0("Unknown_", row_number()),
                  Name),
    orig_longitude = longitude,
    orig_latitude  = latitude
  )

## ---- DBSCAN clustering per farmer (50m) ---------------------------------- ##
Tist_sf <- st_as_sf(
  TistDat6,
  coords = c("longitude", "latitude"),
  crs = 4326
) %>%
  st_transform(32736)  # UTM metres

Tist_sf <- Tist_sf %>%
  group_split(Name) %>%
  map_dfr(function(farmer_df) {
    coords <- st_coordinates(farmer_df)
    farmer_df$grove_cluster <- dbscan(coords, eps = 50, minPts = 1)$cluster
    farmer_df
  })

## ---- Aggregate groves per farmer cluster --------------------------------- ##
Tist_merged <- Tist_sf %>%
  group_by(Name, grove_cluster) %>%
  summarise(
    Name            = first(Name),
    Trees           = sum(Trees, na.rm = TRUE),
    Area_Ha         = sum(Area_Ha, na.rm = TRUE),
    Group_No        = first(Group_No),
    Group_Name      = first(Group_Name),
    Cluster         = first(Cluster),
    Village         = first(Village),
    Group_Center    = first(Group_Center),
    Proj_Area       = first(Proj_Area),
    reg_date        = first(reg_date),
    Admin_Districts = first(Admin_Districts),
    Region          = first(Region),
    orig_longitude  = first(orig_longitude),
    orig_latitude   = first(orig_latitude),
    cluster_size    = n(),  # number of groves in cluster
    geometry        = if (n() == 1) first(geometry) else st_union(geometry),
    .groups         = "drop"
  )

## Convert back to lon/lat: original coords for single-grove clusters,
## mean lon/lat of all points for merged clusters
TistDat7 <- Tist_merged %>%
  st_transform(4326) %>%
  rowwise() %>%
  mutate(
    longitude = if (cluster_size == 1) orig_longitude else mean(st_coordinates(geometry)[, 1]),
    latitude  = if (cluster_size == 1) orig_latitude  else mean(st_coordinates(geometry)[, 2])
  ) %>%
  ungroup() %>%
  select(-geometry, -cluster_size, -orig_longitude, -orig_latitude)

head(TistDat6)
head(TistDat7)


################################################################################
# 4. Pseudonymisation
################################################################################

TistDat7 <- TistDat7 %>%
  mutate(
    Location_ID = paste0("L", dense_rank(paste(latitude, longitude, sep = "_"))),
    Farmer_ID   = paste0("F", dense_rank(Name)),
    Village_ID  = paste0("V", dense_rank(Village)),
    Group_ID    = paste0("G", dense_rank(Group_No)),
    Cluster_ID  = paste0("C", dense_rank(Cluster))
  )

TistDat7_anonymised <- TistDat7 %>%
  select(
    Farmer_ID, Group_ID, Location_ID,
    Trees, reg_date, Area_Ha,
    Village_ID, Cluster_ID,
    latitude, longitude,
    Admin_Districts, Proj_Area, Region
  ) %>%
  distinct()

readr::write_csv(TistDat7_anonymised, "Data/TistDat_cleaned_anonymised.csv")


################################################################################
# 5. LOAD ANONYMISED DATA & DEDUPLICATE
################################################################################

TistUg_Data <- read.csv("Data/TistDat_cleaned_anonymised.csv") %>%
  distinct()

Zero_trees <- TistUg_Data %>%
  filter(Trees == 0)

str(TistUg_Data)

sum(duplicated(TistUg_Data))  # should be 0

## Keep only one row per exact farmer-group-location-grove combination
TistUg_Data <- TistUg_Data %>%
  distinct(Farmer_ID, Group_ID, Location_ID, Trees, Area_Ha, reg_date, .keep_all = TRUE)


################################################################################
# 6. SPATIAL FEATURE CONSTRUCTION
################################################################################

## Add location attributes.
##
## (1) Minimum distance between neighbouring groves — conceptually important
##     for modelling opportunities for peer learning, observation, or
##     informal knowledge spillovers.
##       - Adoption of improved technologies depends on geographic proximity
##         to adopters; learning was strongest from nearby neighbours
##         (Conley & Udry, 2010, American Economic Review, 100(1), 35-69).
##       - Visible outcomes on neighbouring farms encourage adoption; visibility
##         and proximity are crucial for perceived success
##         (Meijer et al., 2015, Agricultural Systems, 135, 83-94).
##       - Adoption is more likely when farmers have geographically proximate
##         peers using the technology
##         (Matuschke & Qaim, 2009, Agricultural Economics, 40(5), 493-505).

## ---- 6.1 Nearest neighbour & first-ring mean distance -------------------- ##

ring_radius <- 2000  # first-ring radius, metres

Location_Coords <- TistUg_Data %>%
  distinct(Location_ID, latitude, longitude) %>%
  rename(Lat = latitude, Lon = longitude)

coords_matrix <- Location_Coords %>%
  select(Lon, Lat) %>%
  as.matrix()
rownames(coords_matrix) <- Location_Coords$Location_ID

dist_matrix <- distm(coords_matrix, fun = distHaversine)  # pairwise Haversine distances (m)
diag(dist_matrix) <- NA  # remove self-distances

min_neighbor_dist <- apply(dist_matrix, 1, min, na.rm = TRUE)

mean_neighbor_first_ring <- apply(dist_matrix, 1, function(x) {
  ring_vals <- x[x <= ring_radius]
  mean(ring_vals, na.rm = TRUE)
})

Neighbor_Distance_Summary <- Location_Coords %>%
  mutate(
    Min_Neighbor_Dist_m   = as.numeric(min_neighbor_dist),
    Mean_FirstRing_Dist_m = as.numeric(mean_neighbor_first_ring)
  ) %>%
  select(Location_ID, Min_Neighbor_Dist_m, Mean_FirstRing_Dist_m)

## ---- 6.2 Distance to nearest forest reserve ------------------------------ ##
##
## Could capture proximity to ecological anchors, microclimatic influencers,
## and/or policy-relevant buffers.
##   - Proximity to forests is often linked to better agroforestry outcomes
##     due to biophysical co-benefits and learning from forest-edge ecologies
##     (Miller et al., 2020, Environmental Evidence, 9(1)).
##   - Forests create cooler, more stable microclimates that can influence
##     surrounding agricultural activity, potentially reducing stress on
##     young seedlings and increasing survival
##     (De Frenne et al., 2019, Nature Ecology & Evolution, 3, 744-749).
##   - Forest loss is lower near protected areas, suggesting greater policy
##     pressure and oversight shaping adjacent agroforestry behaviour
##     (Hansen, Stehman & Potapov, 2010, PNAS, 107(19), 8650-8655).
##   - Proximity to different reserve types correlates with land management
##     and enforcement practices
##     (Nelson & Chomitz, 2011, PLoS ONE, 6(8), e22722).

forest_reserves <- st_read("Data/Uganda_Forest_Reserves/Uganda_Forest_Reserves.shp")

grove_points <- st_as_sf(Location_Coords, coords = c("Lon", "Lat"), crs = 4326)

grove_points_utm    <- st_transform(grove_points, 32636)
forest_reserves_utm <- st_transform(forest_reserves, 32636)

distances <- st_distance(grove_points_utm, forest_reserves_utm)
min_dists <- apply(distances, 1, min)

Location_Dist_Forest <- tibble(
  Location_ID      = Location_Coords$Location_ID,
  Dist_To_Forest_m = as.numeric(min_dists)
)

## ---- 6.3 Join spatial summaries back to the main data --------------------- ##

TistUg_Data <- TistUg_Data %>%
  select(
    -starts_with("Min_Neighbor"),
    -starts_with("Mean_Neighbor"),
    -starts_with("Dist_To_Forest")
  ) %>%
  left_join(Neighbor_Distance_Summary, by = "Location_ID") %>%
  left_join(Location_Dist_Forest, by = "Location_ID") %>%
  mutate(Planted_Tree_Density = Trees / Area_Ha)

write.csv(TistUg_Data, "Data/cleaned_anonymised_TistDat_geodistanced.csv", row.names = FALSE)


################################################################################
# 7. NEIGHBOUR COUNTS BY DISTANCE RING & EXPOSURE INDEX
################################################################################

TistUgGeo_Data <- read.csv("Data/cleaned_anonymised_TistDat_geodistanced.csv") %>%
  distinct()

head(TistUgGeo_Data)

TistUgGeo_Data <- TistUgGeo_Data %>%
  mutate(
    longitude_orig = longitude,
    latitude_orig  = latitude
  )

## Count how many farmers are within 500m, 1000m, ... of every individual farmer
farmers_sf <- st_as_sf(
  TistUgGeo_Data,
  coords = c("longitude", "latitude"),
  crs = 4326
) %>%
  st_transform(32636)  # UTM 36N (Uganda)

buffers <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500)

neighbor_counts <- map_dfc(buffers, function(d) {
  neigh  <- st_is_within_distance(farmers_sf, farmers_sf, dist = d)
  counts <- lengths(neigh) - 1  # exclude self
  tibble(!!paste0("N_within_", d, "m") := counts)
})

farmers_sf <- bind_cols(farmers_sf, neighbor_counts)

## Convert cumulative counts into distance-band counts
farmers_sf <- farmers_sf %>%
  mutate(
    N_0_500m     = N_within_500m,
    N_500_1000m  = N_within_1000m - N_within_500m,
    N_1000_1500m = N_within_1500m - N_within_1000m,
    N_1500_2000m = N_within_2000m - N_within_1500m,
    N_2000_2500m = N_within_2500m - N_within_2000m,
    N_2500_3000m = N_within_3000m - N_within_2500m,
    N_3000_3500m = N_within_3500m - N_within_3000m,
    N_3500_4000m = N_within_4000m - N_within_3500m,
    N_4000_4500m = N_within_4500m - N_within_4000m,
    N_4500_5000m = N_within_5000m - N_within_4500m,
    N_5000_5500m = N_within_5500m - N_within_5000m
  )

summary(farmers_sf$N_within_1000m)
colnames(farmers_sf)

farmers1_sf <- farmers_sf %>%
  dplyr::select(-starts_with("N_within_")) %>%
  st_drop_geometry() %>%
  dplyr::rename(
    latitude  = latitude_orig,
    longitude = longitude_orig
  )

head(farmers1_sf)

## ---- Exposure index -------------------------------------------------------
## Exposure = how many nearby neighbours exist, regardless of how much they
## have planted. Closer neighbours are given greater weight (distance decay).
farmers1_sf <- farmers1_sf %>%
  mutate(
    Total_N =
      N_0_500m + N_500_1000m + N_1000_1500m +
      N_1500_2000m + N_2000_2500m +
      N_2500_3000m + N_3000_3500m +
      N_3500_4000m + N_4000_4500m +
      N_4500_5000m + N_5000_5500m,

    Exposure =
      1.0   * N_0_500m +
      0.8   * N_500_1000m +
      0.6   * N_1000_1500m +
      0.4   * N_1500_2000m +
      0.2   * N_2000_2500m +
      0.1   * N_2500_3000m +
      0.05  * N_3000_3500m +
      0.025 * N_3500_4000m +
      0.012 * N_4000_4500m +
      0.006 * N_4500_5000m +
      0.003 * N_5000_5500m
  )

head(farmers1_sf)

write.csv(farmers1_sf, "Data/cleaned_anonymised_TistDat_geodistanced_neighbors.csv", row.names = FALSE)


################################################################################
# 8. PLANTED TREE DENSITY — SENSITIVITY TESTING (CONSTRUCTED VARIABLE)
################################################################################

TistDat_density <- read.csv("Data/cleaned_anonymised_TistDat_geodistanced_neighbors.csv") %>%
  distinct()

head(TistDat_density)

## Note: planted tree density is sensitive to very small areas, reporting
## errors in area, extreme planters (very high numbers), and nonlinear
## scaling effects.

## ---- 8.1 Filter zero-tree farmers and very small plots -------------------
TistDat_density <- TistDat_density %>%
  filter(Trees > 0, Area_Ha > 0) %>%
  mutate(
    Density_clean = ifelse(is.finite(Planted_Tree_Density), Planted_Tree_Density, NA_real_)
  )

summary(TistDat_density$Trees)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.0    39.0   105.0   324.6   286.0 56896.0

summary(TistDat_density$Area_Ha)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.0010   0.1580   0.3280   0.6899   0.6450 100.0000

summary(TistDat_density$Density_clean)

quantile(
  TistDat_density$Density_clean,
  probs = c(0.5, 0.9, 0.95, 0.97, 0.99, 1),
  na.rm = TRUE
)
#      50%       90%       95%       97%       99%
# 390.4762 1547.4453 2122.7501 2669.1439 5116.8615

## ---- 8.2 Winsorize planted tree density at the 99th percentile -----------
q99 <- quantile(TistDat_density$Density_clean, 0.99, na.rm = TRUE)

TistDat_density_clean <- TistDat_density %>%
  mutate(Density_winsor99 = pmin(Density_clean, q99))

summary(TistDat_density_clean$Density_winsor99)

p_density_winsor99 <- ggplot(TistDat_density_clean, aes(x = Density_winsor99)) +
  geom_histogram(bins = 50, alpha = 0.6) +
  labs(
    x     = "Planted tree density (winsorized at 99th percentile)",
    y     = "Count",
    title = "Distribution of planted tree density"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(size = 10),
    axis.text.x      = element_text(angle = 90, hjust = 1),
    axis.title       = element_text(size = 13)
  )

p_density_winsor99

head(TistDat_density_clean)
nrow(TistDat_density_clean)  ## 33894

ggsave("Output/Manuscript 3 graphs/p_density_winsor99.png",
       p_density_winsor99, dpi = 600, width = 12, height = 6)

## ---- 8.3 Minimum planting-area threshold ----------------------------------
## Small planting-area measurements could be a root cause of density errors;
## check how small plots are and how this could affect density.
quantile(
  TistDat_density_clean$Area_Ha,
  probs = c(0, 0.01, 0.03, 0.05, 0.1),
  na.rm = TRUE
)
#    0%    1%    3%    5%   10%
# 0.001 0.020 0.039 0.052 0.080
## Very small area values are prone to reporting errors and could inflate
## density; setting a minimum planting-area threshold.

TistDat_density_clean <- TistDat_density_clean %>%
  filter(Area_Ha >= 0.03)
## Assuming a horizontal GIS equipment measurement accuracy of +/-10m, a
## conservative minimum mapping area is approximated by the area of a circle
## = 0.0314 ha.

nrow(TistDat_density_clean)  ## 38658

write.csv(TistDat_density_clean,
          "Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity.csv",
          row.names = FALSE)


################################################################################
# 9. NEAR-FAR DOMINANCE
################################################################################

TistDat_NearFar <- read.csv("Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity.csv") %>%
  distinct()

head(TistDat_NearFar)

## Computing Near-Far dominance shares across all distance thresholds.
## Progressive Near-Far framing across cumulative distance bands
## (Near_k = 0 to k*500m, Far_k = k*500m to 5500m).

TistDat_NearFar <- TistDat_NearFar %>%
  mutate(
    ## Near 1: 0-500 m   | Far 1: 500-5500 m
    Near_1 = N_0_500m,
    Far_1  = N_500_1000m + N_1000_1500m + N_1500_2000m + N_2000_2500m +
             N_2500_3000m + N_3000_3500m + N_3500_4000m + N_4000_4500m +
             N_4500_5000m + N_5000_5500m,

    ## Near 2: 0-1000 m  | Far 2: 1000-5500 m
    Near_2 = N_0_500m + N_500_1000m,
    Far_2  = N_1000_1500m + N_1500_2000m + N_2000_2500m + N_2500_3000m +
             N_3000_3500m + N_3500_4000m + N_4000_4500m + N_4500_5000m +
             N_5000_5500m,

    ## Near 3: 0-1500 m  | Far 3: 1500-5500 m
    Near_3 = N_0_500m + N_500_1000m + N_1000_1500m,
    Far_3  = N_1500_2000m + N_2000_2500m + N_2500_3000m + N_3000_3500m +
             N_3500_4000m + N_4000_4500m + N_4500_5000m + N_5000_5500m,

    ## Near 4: 0-2000 m  | Far 4: 2000-5500 m
    Near_4 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m,
    Far_4  = N_2000_2500m + N_2500_3000m + N_3000_3500m + N_3500_4000m +
             N_4000_4500m + N_4500_5000m + N_5000_5500m,

    ## Near 5: 0-2500 m  | Far 5: 2500-5500 m
    Near_5 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m + N_2000_2500m,
    Far_5  = N_2500_3000m + N_3000_3500m + N_3500_4000m + N_4000_4500m +
             N_4500_5000m + N_5000_5500m,

    ## Near 6: 0-3000 m  | Far 6: 3000-5500 m
    Near_6 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m +
             N_2000_2500m + N_2500_3000m,
    Far_6  = N_3000_3500m + N_3500_4000m + N_4000_4500m + N_4500_5000m +
             N_5000_5500m,

    ## Near 7: 0-3500 m  | Far 7: 3500-5500 m
    Near_7 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m +
             N_2000_2500m + N_2500_3000m + N_3000_3500m,
    Far_7  = N_3500_4000m + N_4000_4500m + N_4500_5000m + N_5000_5500m,

    ## Near 8: 0-4000 m  | Far 8: 4000-5500 m
    Near_8 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m +
             N_2000_2500m + N_2500_3000m + N_3000_3500m + N_3500_4000m,
    Far_8  = N_4000_4500m + N_4500_5000m + N_5000_5500m,

    ## Near 9: 0-4500 m  | Far 9: 4500-5500 m
    Near_9 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m +
             N_2000_2500m + N_2500_3000m + N_3000_3500m + N_3500_4000m +
             N_4000_4500m,
    Far_9  = N_4500_5000m + N_5000_5500m,

    ## Near 10: 0-5000 m | Far 10: 5000-5500 m
    Near_10 = N_0_500m + N_500_1000m + N_1000_1500m + N_1500_2000m +
              N_2000_2500m + N_2500_3000m + N_3000_3500m + N_3500_4000m +
              N_4000_4500m + N_4500_5000m,
    Far_10  = N_5000_5500m
  )

head(TistDat_NearFar)

## ---- Dominance shares ------------------------------------------------------
TistDat_NearFar1 <- TistDat_NearFar %>%
  mutate(
    NearFarDominance_1  = ifelse(Near_1  + Far_1  > 0, Near_1  / (Near_1  + Far_1),  NA_real_),
    NearFarDominance_2  = ifelse(Near_2  + Far_2  > 0, Near_2  / (Near_2  + Far_2),  NA_real_),
    NearFarDominance_3  = ifelse(Near_3  + Far_3  > 0, Near_3  / (Near_3  + Far_3),  NA_real_),
    NearFarDominance_4  = ifelse(Near_4  + Far_4  > 0, Near_4  / (Near_4  + Far_4),  NA_real_),
    NearFarDominance_5  = ifelse(Near_5  + Far_5  > 0, Near_5  / (Near_5  + Far_5),  NA_real_),
    NearFarDominance_6  = ifelse(Near_6  + Far_6  > 0, Near_6  / (Near_6  + Far_6),  NA_real_),
    NearFarDominance_7  = ifelse(Near_7  + Far_7  > 0, Near_7  / (Near_7  + Far_7),  NA_real_),
    NearFarDominance_8  = ifelse(Near_8  + Far_8  > 0, Near_8  / (Near_8  + Far_8),  NA_real_),
    NearFarDominance_9  = ifelse(Near_9  + Far_9  > 0, Near_9  / (Near_9  + Far_9),  NA_real_),
    NearFarDominance_10 = ifelse(Near_10 + Far_10 > 0, Near_10 / (Near_10 + Far_10), NA_real_)
  )

head(TistDat_NearFar1)
nrow(TistDat_NearFar1)

## ---- Residual Near-Far dominance (net of Exposure) ------------------------
TistDat_NearFar2 <- TistDat_NearFar1 %>%
  mutate(
    NearFar_resid_1  = resid(lm(NearFarDominance_1  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_2  = resid(lm(NearFarDominance_2  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_3  = resid(lm(NearFarDominance_3  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_4  = resid(lm(NearFarDominance_4  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_5  = resid(lm(NearFarDominance_5  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_6  = resid(lm(NearFarDominance_6  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_7  = resid(lm(NearFarDominance_7  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_8  = resid(lm(NearFarDominance_8  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_9  = resid(lm(NearFarDominance_9  ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude)),
    NearFar_resid_10 = resid(lm(NearFarDominance_10 ~ Exposure, data = TistDat_NearFar1, na.action = na.exclude))
  )

## ---- Sanity checks ---------------------------------------------------------
sum(is.na(TistDat_NearFar2$NearFarDominance_1))
sum(is.na(TistDat_NearFar2$NearFar_resid_1))

cor(
  TistDat_NearFar2$NearFar_resid_1,
  TistDat_NearFar2$Exposure,
  use = "complete.obs"
)  # residuals should be uncorrelated with Exposure

head(TistDat_NearFar2)
nrow(TistDat_NearFar2)


################################################################################
# 10. H4b VARIABLE CONSTRUCTION — NEIGHBOUR-MEAN OUTCOME & DISSIMILARITY
################################################################################
#
# Built for BOTH outcomes used in the main intensity hypotheses (H2/H3/H4):
#   - Density_winsor99  (rate: planting intensity, size-adjusted)
#   - Trees             (raw count: absolute scale of adoption)
#
# Construction-only: this section builds and saves the variables. Model
# integration (H2b/H3b/H4b) happens in a separate script.
#
# Mirrors the existing Near_1..Near_10 cumulative-distance-band logic exactly,
# so NeighborMean_<Outcome>_k lines up one-to-one with Near_k / NearFarDominance_k.
#
# Requires: TistDat_NearFar2 with Location_ID, latitude, longitude,
# Density_winsor99, and Trees already built as above.
################################################################################

## ---- 10.1 Setup ------------------------------------------------------------

OUTCOME_VARS    <- c("Density_winsor99", "Trees")             # outcomes to build similarity variables for
near_thresholds <- seq(500, 5000, by = 500)                    # cumulative distance bands (matches Near_1..Near_10)

farms_sf <- st_as_sf(
  TistDat_NearFar2,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(32636)  # UTM 36N, matching the existing pipeline

## ---- 10.2 Neighbour-mean outcome construction ------------------------------

## For a given outcome and distance, returns the mean of that outcome among
## neighbours within the distance (excluding self).
compute_neighbor_mean <- function(sf_obj, dist_m, outcome) {
  neigh_list <- st_is_within_distance(sf_obj, sf_obj, dist = dist_m)
  vapply(seq_along(neigh_list), function(i) {
    idx <- setdiff(neigh_list[[i]], i)  # exclude self
    if (length(idx) == 0) return(NA_real_)
    mean(outcome[idx], na.rm = TRUE)
  }, numeric(1))
}

## Build NeighborMean_<Outcome>_k for every outcome x threshold combination
for (outcome_name in OUTCOME_VARS) {
  outcome_vec <- farms_sf[[outcome_name]]

  neighbor_means <- map_dfc(seq_along(near_thresholds), function(k) {
    d        <- near_thresholds[k]
    col_name <- paste0("NeighborMean_", outcome_name, "_", k)
    tibble(!!col_name := compute_neighbor_mean(farms_sf, d, outcome_vec))
  })

  TistDat_NearFar2 <- bind_cols(TistDat_NearFar2, neighbor_means)
}

## ---- 10.3 Dissimilarity construction ---------------------------------------
##
## For each outcome x threshold:
##   Dissimilarity_<Outcome>_k    = own - neighbour mean   (signed)
##   AbsDissimilarity_<Outcome>_k = |own - neighbour mean|  (symmetric --
##                                   primary "similarity to neighbours'
##                                   outcomes" measure; low value = high
##                                   similarity)

for (outcome_name in OUTCOME_VARS) {
  own <- TistDat_NearFar2[[outcome_name]]

  for (k in seq_along(near_thresholds)) {
    nbrm_col <- paste0("NeighborMean_", outcome_name, "_", k)
    nbrm     <- TistDat_NearFar2[[nbrm_col]]

    TistDat_NearFar2[[paste0("Dissimilarity_", outcome_name, "_", k)]]    <- own - nbrm
    TistDat_NearFar2[[paste0("AbsDissimilarity_", outcome_name, "_", k)]] <- abs(own - nbrm)
  }
}

## ---- 10.4 Sanity checks -----------------------------------------------------

## Missingness should only come from farms with zero neighbours within a
## given band -- check it's not unexpectedly high, especially at small k.
map_dfr(OUTCOME_VARS, function(outcome_name) {
  map_dfr(seq_along(near_thresholds), function(k) {
    col <- paste0("AbsDissimilarity_", outcome_name, "_", k)
    tibble(
      outcome     = outcome_name,
      k           = k,
      threshold_m = near_thresholds[k],
      n_missing   = sum(is.na(TistDat_NearFar2[[col]])),
      pct_missing = round(100 * mean(is.na(TistDat_NearFar2[[col]])), 2)
    )
  })
})

## Quick distribution check at a mid-range threshold (k = 2, 0-1000m) for
## both outcomes -- confirm nothing looks obviously broken before modelling.
summary(TistDat_NearFar2$AbsDissimilarity_Density_winsor99_2)
summary(TistDat_NearFar2$AbsDissimilarity_Trees_2)

write.csv(TistDat_NearFar2,
          "Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes.csv",
          row.names = FALSE)


################################################################################
# 11. REGISTRATION DURATION
################################################################################

TistDat_regduration <- read.csv("Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes.csv") %>%
  distinct()

head(TistDat_regduration)

ref_year <- 2024

TistDat_regduration <- TistDat_regduration %>%
  dplyr::mutate(Years_since_reg = ref_year - reg_date)

head(TistDat_regduration)


################################################################################
# 12. ADD SUBCOUNTY (SPATIAL JOIN)
################################################################################

subcounty <- st_read("C:/Users/ae474/OneDrive - University of Exeter/PhD work in progress/Data/R Large Data/UGshapefiles/uga_admbnda_adm4_ubos_20200824.shp")

Tistfarmers_sf <- st_as_sf(
  TistDat_regduration,
  coords = c("longitude", "latitude"),
  crs = 4326
)

## Confirm CRS matches
st_crs(subcounty)
st_crs(Tistfarmers_sf)

names(subcounty)  ## subcounty names = ADM4_EN

subcounty <- st_make_valid(subcounty)  ## repair subcounty geometry

farmers_with_subcounty <- st_join(
  Tistfarmers_sf,
  subcounty[, "ADM4_EN"],  # change if column name differs
  join = st_within
) %>%
  rename(Subcounty = ADM4_EN) %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude  = st_coordinates(.)[, 2]
  )

TistDat_regduration <- farmers_with_subcounty %>%
  st_drop_geometry()

table(is.na(TistDat_regduration$Subcounty))  ## check for NAs

head(TistDat_regduration)

write.csv(TistDat_regduration,
          "Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv",
          row.names = FALSE)


################################################################################
# 13. JOIN NPHC SUBCOUNTY POPULATION DATA (BUSHENYI & SOROTI)
################################################################################

NPHC_subcounty <- read.csv("Data/BushSoroti_subcountyNPHC.csv") %>%
  distinct()

BushSoroti_TistDat_regduration <- read.csv("Data/cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv") %>%
  distinct() %>%
  filter(Proj_Area %in% c("Bushenyi", "Soroti"))

head(NPHC_subcounty)
head(BushSoroti_TistDat_regduration)

B_NPHC_subcounty <- NPHC_subcounty %>%
  filter(Admin_Districts %in% c(
    "SHEEMA", "MITOOMA", "RUBIRIZI", "BUSHENYI",
    "BUHWEJU", "KITAGWENDA", "IBANDA", "NTUNGAMO",
    "KAMWENGE", "KASESE", "RWAMPARA", "KYENJOJO"
  ))
head(B_NPHC_subcounty)
n_distinct(B_NPHC_subcounty$Subcounty)

S_NPHC_subcounty <- NPHC_subcounty %>%
  filter(Admin_Districts %in% c(
    "SOROTI", "ALEBTONG", "KALAKI",
    "AMURIA", "SERERE", "KAPELEBYONG"
  ))
n_distinct(S_NPHC_subcounty$Subcounty)

## Join NPHC_subcounty to TistDat_regduration on exact name matches only.
## Standardise case first.
NPHC_subcounty_clean <- NPHC_subcounty %>%
  mutate(Subcounty = toupper(Subcounty)) %>%
  mutate(Total_Subcounty_popn. = as.numeric(gsub(",", "", Total_Subcounty_popn.)))

BushSoroti_TistDat_regduration_clean <- BushSoroti_TistDat_regduration %>%
  mutate(Subcounty = toupper(Subcounty))

Tist_joined <- BushSoroti_TistDat_regduration_clean %>%
  left_join(NPHC_subcounty_clean, by = "Subcounty")

head(Tist_joined)

write.csv(Tist_joined,
          "Data/BushSoroti_cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv",
          row.names = FALSE)


################################################################################
# 14. JOIN NPHC HOUSEHOLD COMMUNICATION DATA
################################################################################

BushSoroti_TISTsubcountypopn <- read.csv("Data/BushSoroti_cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered.csv") %>%
  distinct()

NPHC_communication <- read.csv("Data/NPHC2024_household_communication.csv") %>%
  distinct()

head(NPHC_communication)

NPHC_communication <- NPHC_communication %>%
  mutate(
    Total_Households    = as.numeric(gsub(",", "", Total_Households)),
    Radio                = as.numeric(gsub(",", "", Radio)),
    Word_of_Mouth        = as.numeric(gsub(",", "", Word_of_Mouth)),
    Phone_calls          = as.numeric(gsub(",", "", Phone_calls)),
    Community_meetings   = as.numeric(gsub(",", "", Community_meetings)),
    Community_Announcer  = as.numeric(gsub(",", "", Community_Announcer))
  )

Tist_commuication_joined <- BushSoroti_TISTsubcountypopn %>%
  left_join(NPHC_communication, by = "Subcounty")

head(Tist_commuication_joined)

write.csv(Tist_commuication_joined,
          "Data/BushSoroti_cleaned_anonymised_TistDat_geodistanced_neighbors_winsoriseddensity_NearFarRes_InfluenceFieldRes_dateregistered_communicationoptions.csv",
          row.names = FALSE)


################################################################################
# 15. DIAGNOSTIC — MISSINGNESS IN NeighborMean / Dissimilarity VARIABLES
################################################################################
#
# Per threshold (k = 1..10), per outcome (Density_winsor99, Trees), by Proj_Area.
#
# Purpose: check whether missingness (farms with zero neighbours within a
# given band) is (a) small enough to ignore, (b) concentrated at small k only,
# and (c) evenly distributed across Bushenyi vs Soroti -- or whether it's
# disproportionately dropping one site, which would bias threshold comparisons.
#
# Run on TistDat_NearFar2 (or Tist_commuication_joined / whatever object
# currently holds the NeighborMean_/Dissimilarity_ columns).
################################################################################

OUTCOME_VARS    <- c("Density_winsor99", "Trees")
near_thresholds <- seq(500, 5000, by = 500)

## ---- 15.1 Overall missingness per outcome x threshold ---------------------

missing_overall <- map_dfr(OUTCOME_VARS, function(outcome_name) {
  map_dfr(seq_along(near_thresholds), function(k) {
    col <- paste0("AbsDissimilarity_", outcome_name, "_", k)
    tibble(
      outcome     = outcome_name,
      k           = k,
      threshold_m = near_thresholds[k],
      n_total     = nrow(TistDat_NearFar2),
      n_missing   = sum(is.na(TistDat_NearFar2[[col]])),
      pct_missing = round(100 * mean(is.na(TistDat_NearFar2[[col]])), 2)
    )
  })
})

cat("\n=== Overall missingness by threshold and outcome ===\n")
print(missing_overall, n = 30)

## ---- 15.2 Missingness split by Proj_Area (key check for site-level bias) --

missing_by_area <- map_dfr(OUTCOME_VARS, function(outcome_name) {
  map_dfr(seq_along(near_thresholds), function(k) {
    col <- paste0("AbsDissimilarity_", outcome_name, "_", k)
    TistDat_NearFar2 %>%
      group_by(Proj_Area) %>%
      summarise(
        outcome     = outcome_name,
        k           = k,
        threshold_m = near_thresholds[k],
        n_total     = n(),
        n_missing   = sum(is.na(.data[[col]])),
        pct_missing = round(100 * mean(is.na(.data[[col]])), 2),
        .groups = "drop"
      )
  })
})

cat("\n=== Missingness by threshold, outcome, and Proj_Area ===\n")
print(missing_by_area, n = 100)

## ---- 15.3 Plot: missingness curve across thresholds, split by area --------
##      (one panel per outcome, lines coloured by Proj_Area)

missing_plot <- ggplot(
  missing_by_area,
  aes(x = threshold_m, y = pct_missing, color = Proj_Area)
) +
  geom_line() +
  geom_point() +
  facet_wrap(~ outcome) +
  labs(
    title    = "Missingness in AbsDissimilarity by threshold and project area",
    subtitle = "Missing = farm has zero neighbours within that cumulative distance band",
    x        = "Distance threshold (m)",
    y        = "% missing",
    color    = "Project area"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.title       = element_text(size = 13)
  )

missing_plot

## ---- 15.4 Formal check: is missingness independent of Proj_Area? ----------
##
## Chi-square test per threshold (using Density_winsor99 as the reference
## outcome, since Density_winsor99/Trees share the same underlying spatial
## missingness pattern -- missing is about neighbour COUNT, not the outcome
## value itself, so one test per threshold is sufficient).

chisq_by_threshold <- map_dfr(seq_along(near_thresholds), function(k) {
  col <- paste0("AbsDissimilarity_Density_winsor99_", k)
  tab <- table(
    TistDat_NearFar2$Proj_Area,
    is.na(TistDat_NearFar2[[col]])
  )
  test <- suppressWarnings(chisq.test(tab))
  tibble(
    k           = k,
    threshold_m = near_thresholds[k],
    chisq_stat  = unname(test$statistic),
    p_value     = test$p.value
  )
})

cat("\n=== Chi-square test: missingness independent of Proj_Area? ===\n")
cat("(low p-value = missingness DOES differ significantly by area at that threshold)\n")
print(chisq_by_threshold, n = 10)

## ---- 15.5 Isolation summary -------------------------------------------------
##
## How many farms are missing at EVERY threshold (i.e. isolated even at
## 5000m) vs. only missing at small k (i.e. just locally sparse)?

isolation_summary <- TistDat_NearFar2 %>%
  mutate(
    missing_at_1  = is.na(AbsDissimilarity_Density_winsor99_1),
    missing_at_10 = is.na(AbsDissimilarity_Density_winsor99_10)
  ) %>%
  summarise(
    n_total                   = n(),
    n_missing_500m_only       = sum(missing_at_1 & !missing_at_10),
    n_missing_even_at_5000m   = sum(missing_at_10),
    pct_missing_even_at_5000m = round(100 * mean(missing_at_10), 2)
  )

cat("\n=== Isolation summary ===\n")
print(isolation_summary)

isolation_by_area <- TistDat_NearFar2 %>%
  mutate(missing_at_10 = is.na(AbsDissimilarity_Density_winsor99_10)) %>%
  group_by(Proj_Area) %>%
  summarise(
    n_total                   = n(),
    n_missing_even_at_5000m   = sum(missing_at_10),
    pct_missing_even_at_5000m = round(100 * mean(missing_at_10), 2),
    .groups = "drop"
  )

cat("\n=== Isolation (missing even at 5000m) by Proj_Area ===\n")
print(isolation_by_area)

## ---- 15.6 Save diagnostics --------------------------------------------------

ggsave("Output/Manuscript 3 graphs/H4b_missingness_diagnostic.png",
       missing_plot, dpi = 600, width = 10, height = 6)

write.csv(missing_by_area, "Output/H4b_missingness_by_area_threshold.csv", row.names = FALSE)
write.csv(chisq_by_threshold, "Output/H4b_missingness_chisq_by_threshold.csv", row.names = FALSE)

################################################################################
# END OF SCRIPT
################################################################################
