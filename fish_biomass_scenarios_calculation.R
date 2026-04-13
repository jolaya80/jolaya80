# This Code does:
# Load coral scenario rasters, prepare depth and fishing-zone rasters, then join empirical fish survey points to habitat + zone and summarise biomass by habitat and zone.
#

# Julian Olaya-Restrepo 2026


# 0. Setup: libraries and directories
# - Load packages used across the script (spatial, data wrangling, plotting, IO).
#   Keep startup messages quiet for cleaner logs.
suppressPackageStartupMessages({
  library(terra)    # raster-like spatial data (fast, modern)
  library(sf)       # vector spatial data (simple features)
  library(dplyr)    # data wrangling (tibble verbs)
  library(tidyr)    # reshaping data (wide <-> long)
  library(ggplot2)  # plotting
  library(writexl)  # write simple Excel files
  library(readr)    # fast CSV read/write
})

# - Ensure strings are not converted to factors by default.
options(stringsAsFactors = FALSE)

# Base project directory (change to your machine/project as needed)
proj_dir  <- "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models"

# Sub-directories for each model component
fish_dir  <- file.path(proj_dir, "Fish_biomass")
coral_dir <- file.path(proj_dir, "coral_cover_modeling")   # adjust if needed

# Output directories inside fish project
out_dir           <- file.path(fish_dir, "outputs")
tabplot_out_dir   <- file.path(out_dir, "tables_plots")     
spatial_out_dir   <- file.path(out_dir, "spatial_files")

# Create directories if they don't exist. recursive = TRUE makes parent folders as needed.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabplot_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(spatial_out_dir, showWarnings = FALSE, recursive = TRUE)

# Record session information (packages and versions) to the outputs folder for reproducibility.
# This is important for biology / fisheries modeling so others can reproduce analyses later.
session_file <- file.path(out_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = session_file)


############################################################
# 1. Load coral scenario rasters and depth
############################################################

# # - Load coral scenario rasters (baseline + scenarios), load bathymetry,
# #  align bathymetry to coral raster grid, and save the aligned depth raster.

# Paths (use variables so it's easy to change)
coral_raster_dir <- file.path(coral_dir, "06_outputs", "two_scenarios")

# Check files exist before trying to read them (fail early with helpful message)
required_coral_files <- c("baseline_class.tif", 
                          "phase1_class.tif",
                          "phase1_bleach_class.tif", 
                          "phase2_class.tif",
                          "bleaching_only_class.tif")
missing_files <- required_coral_files[!file.exists(file.path(coral_raster_dir, required_coral_files))]
if (length(missing_files) > 0) {
  stop("Missing coral raster files: ", paste(missing_files, collapse = ", "),
       ". Please run the coral workflow or correct paths.")
}

# baseline_class.tif: 1 = Degraded (<15%), 2 = Healthy (>=15%)
# phase1_class.tif, bleach_class.tif, phase2_class.tif

# Read rasters (terra::rast is fast and preserves CRS/resolution metadata)
baseline_class      <- rast(file.path(coral_raster_dir, "baseline_class.tif"))
phase1_class        <- rast(file.path(coral_raster_dir, "phase1_class.tif"))
bleach_class        <- rast(file.path(coral_raster_dir, "phase1_bleach_class.tif"))
phase2_class        <- rast(file.path(coral_raster_dir, "phase2_class.tif"))
bleaching_only_class <- rast(file.path(coral_raster_dir, "bleaching_only_class.tif"))


###############
# Depth raster: ensure the file exists and then align to coral grid

# use this line only if the file was previously created in previous sesions and avoid lines 86 - 118
# depth_res <- rast("C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/spatial_files/depth_res_to_coral.tif")

depth_path <- "G:/Shared drives/NSF CoPE internal/GIS_CoPE/GIS_Belize/2_model_inputs_belize/coral_reef_modeling/05_Preparation_spatial_predictors/10m/depth.tif"
if (!file.exists(depth_path)) {
  stop("Depth raster not found at: ", depth_path)
}
depth_raw <- rast(depth_path)

#Check CRS and reproject depth raster if needed (terra)
crs_depth <- crs(depth_raw)
crs_coral  <- crs(baseline_class)

if (is.na(crs_depth) || crs_depth == "") {
  stop("Depth raster has no CRS. Set it before reprojecting (e.g. crs(depth_raw) <- 'EPSG:4326' or correct EPSG).")
}
if (is.na(crs_coral) || crs_coral == "") {
  stop("Baseline coral raster has no CRS. Inspect source raster.")
}

# If CRS differ, reproject depth to match coral raster (also aligns to template grid)
if (!identical(crs_depth, crs_coral)) {
  message("Reprojecting and aligning depth raster to coral raster CRS/grid...")
  depth_raw <- project(depth_raw, baseline_class, method = "bilinear")
}


# Crop and resample depth to coral grid
depth_crop <- crop(depth_raw, baseline_class)
depth_res  <- resample(depth_crop, baseline_class, method = "bilinear")

# 3) Ensure output directory exists, then save aligned depth raster
dir.create(spatial_out_dir, showWarnings = FALSE, recursive = TRUE)
writeRaster(depth_res,
            filename = file.path(spatial_out_dir, "depth_res_to_coral.tif"),
            overwrite = TRUE)


# 4) Quick checks: print target coral resolution and a small visual check if desired
print(res(baseline_class))
# optional quick visual check
# plot(baseline_class, main = "Baseline coral raster")
# plot(depth_res, add = TRUE)


# Quick checks
message("Coral raster resolution (x, y): ", paste(res(baseline_class), collapse = ", "))
print(baseline_class)
print(is.lonlat(baseline_class))   # FALSE = projected (meters)
print(freq(baseline_class, digits = 6))   # class counts
print(depth_res)
print(global(depth_res, fun = c("min","max","mean"), na.rm = TRUE))

############################################################
# 2. Load fishing zones and rasterize to zone_id
############################################################

# 2. Load fishing zones and rasterize to zone_id
zones_shp <- "G:/Shared drives/NSF CoPE internal/GIS_CoPE/GIS_Belize/_source_data_belize/fishing_zones/Belize Managed Access Areas shape files/Managed Access Areas shape files/ma2016.shp"
if (!file.exists(zones_shp)) stop("Fishing zones shapefile not found: ", zones_shp)

fishing_zones <- st_read(zones_shp, quiet = TRUE)

message("Polygons read: ", nrow(fishing_zones))
message("First few names:")
print(head(fishing_zones$Name))

# Reproject to coral CRS if needed
if (st_crs(fishing_zones) != st_crs(baseline_class)) {
  fishing_zones_utm <- st_transform(fishing_zones, crs(baseline_class))
} else {
  fishing_zones_utm <- fishing_zones
}

# Create numeric zone_id from Name in a reproducible way
fishing_zones_utm$zone_id <- as.numeric(factor(fishing_zones_utm$Name))

# Save reprojected polygons for record
out_gpkg <- file.path(fish_dir, "01_shp", "fishing_zones_utm.gpkg")
dir.create(dirname(out_gpkg), showWarnings = FALSE, recursive = TRUE)
st_write(fishing_zones_utm, out_gpkg, layer = "fishing_zones", delete_layer = TRUE, quiet = TRUE)


# Convert to terra vector and rasterize to coral grid
zones_vect <- vect(fishing_zones_utm)
zones_raster <- rasterize(zones_vect, baseline_class, field = "zone_id", background = NA)
names(zones_raster) <- "zone_id"

# Checks after rasterize
message("Unique zone IDs in raster (excluding NA):")
print(sort(na.omit(unique(values(zones_raster)))))


# Make a lookup table: numeric zone_id -> Name
zone_table <- fishing_zones_utm %>%
  st_drop_geometry() %>%
  dplyr::select(zone_id, Zone = Name) %>%
  distinct()

# Quick check
print(unique(values(zones_raster)))
print(zone_table)

############################################################
# 3. Load empirica fish biomass and build biomass_by_habitat_zone
#    (tCommBioma in g/100 m^2)
############################################################

fish_points <- st_read(
  file.path(fish_dir, "01_shp/3_fish_bz_data.shp")
)

# Ensure fish points are in same CRS as coral rasters
fish_points <- st_transform(fish_points, crs(baseline_class))

# Start habitat classification table
# Only keep key fields for biomass
fish_habitat <- fish_points %>%
  dplyr::select(Name, Latitude, Longitude, tCommBioma, tCORALavg_, geometry) %>%
  mutate(
    tCommBioma = as.numeric(tCommBioma),
    tCORALavg_ = as.numeric(tCORALavg_),
    Habitat = NA_character_
  )

# 3.1 Assign coral reef condition from monitoring tCORALavg_ first (>=15% = Healthy)
fish_habitat <- fish_habitat %>%
  mutate(
    Habitat = if_else(!is.na(tCORALavg_) & tCORALavg_ >= 15,
                      "Coral Reef Healthy",
                      if_else(!is.na(tCORALavg_), "Coral Reef Degraded", NA_character_))
  )

# 3.2 Fallback to coral baseline raster where Habitat still NA
fish_vect_for_extract <- vect(fish_habitat)
coral_vals <- terra::extract(baseline_class, fish_vect_for_extract)[, 2]  # second column is raster values
fish_habitat$Coral_baseline <- coral_vals  # 1 or 2

fish_habitat <- fish_habitat %>%
  mutate(
    Habitat = ifelse(is.na(Habitat) & !is.na(Coral_baseline),
                     ifelse(Coral_baseline == 2, "Coral Reef Healthy", "Coral Reef Degraded"),
                     Habitat)
  )

# QC counts for how habitat was assigned
qc_counts <- fish_habitat %>%
  mutate(
    source = case_when(
      !is.na(tCORALavg_) & !is.na(Coral_baseline) ~ "Monitoreo + Raster",
      !is.na(tCORALavg_) &  is.na(Coral_baseline) ~ "Solo Monitoreo",
      is.na(tCORALavg_) & !is.na(Coral_baseline) ~ "Solo Raster",
      TRUE                                       ~ "Ninguna Fuente"
    )
  ) %>%
  st_drop_geometry() %>%
  count(source)
message("Habitat assignment sources (counts):")
print(qc_counts)

######
# 3.3 Assign fishing zones to fish points (by polygon overlay)
fish_habitat <- st_join(fish_habitat, fishing_zones_utm %>% dplyr::select(zone_id, Name), left = TRUE)

# Rename for clarity
if ("Name.x" %in% names(fish_habitat)) fish_habitat <- rename(fish_habitat, Site = Name.x)
if ("Name.y" %in% names(fish_habitat)) fish_habitat <- rename(fish_habitat, Zone = Name.y)
# If join used only Name, it may already be "Name"
if (!"Zone" %in% names(fish_habitat) & "Name" %in% names(fish_habitat)) {
  fish_habitat <- rename(fish_habitat, Zone = Name)
}

# 3.4 Summarize biomass by Habitat and Zone
# tCommBioma is g/100 m² at each survey point
biomass_by_habitat_zone <- fish_habitat %>%
  st_drop_geometry() %>%
  group_by(Zone, Habitat) %>%
  summarise(
    mean_biomass_g100m2 = mean(tCommBioma, na.rm = TRUE),
    n = sum(!is.na(tCommBioma)),
    .groups = "drop"
  )

# Save this table (important for methods)
write_xlsx(biomass_by_habitat_zone, path = file.path(tabplot_out_dir, "biomass_by_habitat_zone.xlsx"))
message("Saved biomass_by_habitat_zone to: ", file.path(tabplot_out_dir, "biomass_by_habitat_zone.xlsx"))

library(readxl)
# upload table if starting from a new sesion
biomass_by_habitat_zone <- read_xlsx(path = file.path(tabplot_out_dir, "biomass_by_habitat_zone.xlsx"))

## box plot for visualization
library(ggplot2)
library(dplyr)

# Prepare plot data and boxplot (unchanged structure, tidy pipes)
fish_plot <- fish_habitat %>%
  filter(!is.na(Habitat), !is.na(Zone)) %>%
  mutate(Habitat = factor(Habitat, levels = c("Coral Reef Degraded", "Coral Reef Healthy")))


mean_biomass <- fish_habitat %>%
  filter(!is.na(Habitat), !is.na(Zone)) %>%
  group_by(Zone, Habitat) %>%
  summarise(mean_biomass_g100m2 = mean(tCommBioma, na.rm = TRUE), .groups = "drop")


p_box <- ggplot(fish_plot, aes(x = Zone, y = tCommBioma, fill = Habitat)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Habitat),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, alpha = 0.8) +
  geom_point(data = mean_biomass, aes(x = Zone, y = mean_biomass_g100m2, group = Habitat),
             position = position_dodge(width = 0.8),
             shape = 23, size = 3, fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = c("Coral Reef Degraded" = "#E60000", "Coral Reef Healthy"  = "#07F5F8")) +
  scale_color_manual(values = c("Coral Reef Degraded" = "#E60000", "Coral Reef Healthy"  = "#07F5F8")) +
  labs(x = NULL, y = "Fish biomass (g / 100 m^2)", fill = "Reef condition", color = "Reef condition") +
  theme_classic(base_size = 14) +
  theme(text = element_text(family = "sans", size = 12),
        axis.title = element_text(size = 11),
        axis.text  = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))

p_box

# Export boxplot (Frontiers size approximation)
one_col_in  <- 85 / 25.4   # 1 column in inches
two_col_in  <- 170 / 25.4  # 2 columns in inches
ggsave(filename = file.path(out_dir, "biomass_fishing_zone_boxplot.tiff"),
       plot = p_box, device = "tiff",
       width = two_col_in, height = two_col_in * 0.6, units = "in",
       dpi = 300, compression = "lzw")

############################################################
# 4. Build biomass lookup: coral degraded vs healthy per zone
############################################################

biomass_lookup <- biomass_by_habitat_zone %>%
  filter(Habitat %in% c("Coral Reef Degraded", "Coral Reef Healthy")) %>%
  mutate(
    class   = ifelse(Habitat == "Coral Reef Degraded", 1, 2),   # 1=degraded, 2=healthy
    biomass = mean_biomass_g100m2                                # g/100m2
  ) %>%
  group_by(Zone, class) %>%
  summarise(
    biomass = mean(biomass, na.rm = TRUE),
    .groups = "drop"
  )
print(biomass_lookup)

all_zones <- na.omit(unique(biomass_lookup$Zone))

############################################################
# 5. Helper functions and constants
############################################################

# Pixel area in m^2 (explicit)
res_xy <- res(baseline_class)
pixel_area_m2 <- prod(res_xy)  # should match user's raster resolution
message("Pixel area (m^2): ", pixel_area_m2)

# Convert g/100 m² to g per pixel: multiply by (pixel_area / 100)
g100m2_to_g_per_pixel <- function(x) x * (pixel_area_m2 / 100)

# Convert g to kg and tons for summaries
g_to_kg <- function(x) x / 1000
kg_to_ton <- function(x) x / 1000

# Wrapper for zonal stats
safe_zonal <- function(raster_layer, zones, name, fun_type = "sum") {
  df <- zonal(raster_layer, zones, fun = fun_type, na.rm = TRUE) %>%
    as.data.frame()
  names(df) <- c("zone_id", name)
  df
}

############################################################
# 6. Baseline: per-pixel fish biomass raster (g/100m2)
############################################################

# Load the raster file back if starting from new session
# baseline_biomass_pixel <- rast(file.path(spatial_out_dir, "baseline_fish_biomass_g_per_pixel.tif"))

# Initialize raster with NA
baseline_biomass_pixel <- baseline_class * NA
names(baseline_biomass_pixel) <- "baseline_biomass_g_per_pixel"

for (z in all_zones) {
  if (is.na(z)) next
  
  zone_row <- zone_table %>% filter(Zone == z)
  if (nrow(zone_row) == 0) next
  zid <- zone_row$zone_id
  
  zone_pixels <- zones_raster == zid
  deg_pixels     <- zone_pixels & (baseline_class == 1)
  healthy_pixels <- zone_pixels & (baseline_class == 2)
  
  bm_deg <- biomass_lookup %>% filter(Zone == z, class == 1) %>% pull(biomass)
  bm_healthy <- biomass_lookup %>% filter(Zone == z, class == 2) %>% pull(biomass)
  
  if (length(bm_deg) == 0)     bm_deg <- NA
  if (length(bm_healthy) == 0) bm_healthy <- NA
  
  # Convert mean (g / 100 m^2) -> g per pixel
  bm_deg_g_per_pixel     <- if (!is.na(bm_deg))     g100m2_to_g_per_pixel(bm_deg)     else NA_real_
  bm_healthy_g_per_pixel <- if (!is.na(bm_healthy)) g100m2_to_g_per_pixel(bm_healthy) else NA_real_
  
  baseline_biomass_pixel[deg_pixels]     <- bm_deg_g_per_pixel
  baseline_biomass_pixel[healthy_pixels] <- bm_healthy_g_per_pixel
}

writeRaster(baseline_biomass_pixel,
            filename = file.path(spatial_out_dir, "baseline_fish_biomass_g_per_pixel.tif"),
            overwrite = TRUE)

############################################################
# 7. Phase 1: restoration (4.9× on improved pixels)
############################################################

# Load the raster file back if starting from new session
# phase1_biomass_pixel <- rast(file.path(spatial_out_dir, "phase1_fish_biomass_g100m2.tif"))


phase1_biomass_pixel <- baseline_class * NA
names(phase1_biomass_pixel) <- "phase1_biomass_g_per_pixel"

for (z in all_zones) { # Iterates over each element z in the vector all_zones
  if (is.na(z)) next   # If the current zone value is NA, skip
  
  zone_row <- zone_table %>% filter(Zone == z) # From zone_table, select the row corresponding to the current zone z
  if (nrow(zone_row) == 0) next # If no matching row exists for this zone, skip
  zid <- zone_row$zone_id # Pulls the numeric raster ID used to identify this zone in zones_raster
  
  zone_pixels <- zones_raster == zid # Creates a logical raster (TRUE where the raster pixel belongs to zone z, FALSE elsewhere)
  
  deg_base     <- zone_pixels & (baseline_class == 1) # pixels degraded at baseline
  healthy_base <- zone_pixels & (baseline_class == 2) # pixels healthy at baseline
  
  deg_phase1     <- zone_pixels & (phase1_class == 1) # degraded in Phase 1
  healthy_phase1 <- zone_pixels & (phase1_class == 2) # healthy in Phase 1
  
  improved_pixels <- deg_base & healthy_phase1 # Identify improved pixels (captures transition / recovery)
  
  bm_deg <- biomass_lookup %>% # Look up biomass values for degraded class
    filter(Zone == z, class == 1) %>% 
    pull(biomass)
  
  bm_healthy <- biomass_lookup %>% # Look up biomass values for healthy class
    filter(Zone == z, class == 2) %>% 
    pull(biomass)
  
  # Handle missing biomass entries
  if (length(bm_deg) == 0)     bm_deg <- NA
  if (length(bm_healthy) == 0) bm_healthy <- NA
  
  # Convert biomass units to per-pixel grams
  # Convert from grams per 100 m² to grams per pixel
  bm_deg_g_per_pixel     <- if (!is.na(bm_deg))     g100m2_to_g_per_pixel(bm_deg)     else NA_real_
  bm_healthy_g_per_pixel <- if (!is.na(bm_healthy)) g100m2_to_g_per_pixel(bm_healthy) else NA_real_
  
  # 1) Pixels still degraded (baseline & Phase 1)
  # For pixels that:
  # Were degraded at baseline and
  # Are still degraded in Phase 1
  # Assign degraded biomass per pixel.
  phase1_biomass_pixel[deg_phase1 & deg_base] <- bm_deg_g_per_pixel
  
  # 2) Pixels still healthy
  # For pixels that:
  # Were healthy at baseline and
  # Are still healthy in Phase 1
  # Assign healthy biomass per pixel
  phase1_biomass_pixel[healthy_phase1 & healthy_base] <- bm_healthy_g_per_pixel
  
  # 3) Improved pixels: degraded -> healthy, 4.9× from degraded level (applied to per-pixel grams)
  # For pixels that transitioned
  # Degraded → Healthy
  # Assign biomass as:
  # 4.9 × degraded biomass per pixel
  phase1_biomass_pixel[improved_pixels] <- bm_deg_g_per_pixel * 4.9
}

writeRaster(phase1_biomass_pixel,
            filename = file.path(spatial_out_dir, "phase1_fish_biomass_g_per_pixel.tif"),
            overwrite = TRUE)

# phase1_biomass_pixel <- rast("C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/spatial_files/phase1_fish_biomass_g100m2.tif")

############################################################
# 8. Bleaching: apply 57% decline (SF = 0.43)
############################################################

phase1_bleach_biomass_pixel <- phase1_biomass_pixel

# Pixels that bleached: Healthy in Phase 1 -> Degraded after bleaching
bleach_loss_pixels <- (phase1_class == 2) & (bleach_class == 1)

# Apply SF = 0.43 (retain 43% of biomass)
phase1_bleach_biomass_pixel[bleach_loss_pixels] <-
  phase1_biomass_pixel[bleach_loss_pixels] * 0.43

writeRaster(phase1_bleach_biomass_pixel,
            filename = file.path(spatial_out_dir, "phase1_bleach_fish_biomass_g_per_pixel.tif"),
            overwrite = TRUE)


## ELIMINAR
phase1_bleach_biomass_pixel <- rast("C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/spatial_files/phase1_bleach_fish_biomass_g_per_pixel.tif")

############################################################
# 8b. Bleaching Only: apply 57% decline to baseline (SF = 0.43)
#     Counterfactual: bleaching applied directly to baseline
#     (no Phase 1 restoration). Healthy → Degraded pixels lose
#     57 % of their baseline biomass; all other pixels retain
#     their baseline values.
############################################################

bleaching_only_biomass_pixel <- baseline_biomass_pixel
names(bleaching_only_biomass_pixel) <- "bleaching_only_biomass_g_per_pixel"

# Pixels that were Healthy in baseline and became Degraded in bleaching_only
# 1 = Degraded (<15%)
# 2 = Healthy (>=15%)
bleaching_only_loss_pixels <- (baseline_class == 2) & (bleaching_only_class == 1)

# Apply SF = 0.43 (retain 43% of biomass on those pixels)
bleaching_only_biomass_pixel[bleaching_only_loss_pixels] <-
  baseline_biomass_pixel[bleaching_only_loss_pixels] * 0.43

writeRaster(bleaching_only_biomass_pixel,
            filename = file.path(spatial_out_dir, "bleaching_only_fish_biomass_g_per_pixel.tif"),
            overwrite = TRUE)

# ELIMINAR
bleaching_only_biomass_pixel <- rast("C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/spatial_files/bleaching_only_fish_biomass_g_per_pixel.tif")

############################################################
# 8c. Bleaching Only – Phase 2 (~5 Years): partial fish recovery
#     No coral restoration, but deep-reef fish biomass recovers +12%
#     on pixels that bleached (Healthy baseline → Degraded in bleaching_only).
#     Shallow bleached pixels and never-bleached pixels are unchanged.
############################################################

# Initialize from the post-bleach Bleaching Only raster
bleaching_only_phase2_biomass_pixel <- bleaching_only_biomass_pixel
names(bleaching_only_phase2_biomass_pixel) <- "bleaching_only_phase2_biomass_g_per_pixel"

# Depth mask: pixels deeper than 10 m (depth values are negative; <= -10 means >= 10 m deep)
deep_mask <- depth_res <= -10

# Deep recovery: bleached pixels (baseline Healthy → bleaching_only Degraded) that are deep
bleaching_only_deep_recovery <- bleaching_only_loss_pixels & deep_mask

# Apply +12% to those deep pixels relative to the post-bleach biomass
bleaching_only_phase2_biomass_pixel[bleaching_only_deep_recovery] <-
  bleaching_only_biomass_pixel[bleaching_only_deep_recovery] * 1.12

# Note: shallow bleached pixels keep their bleaching_only_biomass_pixel values;
#       never-bleached pixels retain baseline values (inherited from initialisation).

writeRaster(bleaching_only_phase2_biomass_pixel,
            filename = file.path(spatial_out_dir, "bleaching_only_phase2_fish_biomass_g_per_pixel.tif"),
            overwrite = TRUE)
# ELIMINAR
bleaching_only_phase2_biomass_pixel <- rast("C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/spatial_files/bleaching_only_phase2_fish_biomass_g_per_pixel.tif")

############################################################
# 9. Phase 2: continuation of restoration (6×)
#    + depth-based partial recovery (+12%) for deep bleached reefs
############################################################
# Phase 2 biomass assignment (Methods-compliant)
# - Start from post-bleach biomass, then:
#   * For baseline-degraded pixels that are healthy in Phase 2 (either recovered during Phase2
#     or restored in Phase1 and remain healthy), apply scaling factor S2 = 6.9.
#   * Pixels that are healthy across all scenarios retain baseline biomass values.
#   * Pixels degraded by bleaching and remaining degraded after Phase 2:
#       - deep (>10 m): +12% relative to post-bleach biomass
#       - shallow (<10 m): keep post-bleach biomass

# 1) Initialize phase2 raster from the post-bleach raster (keep values by default)
phase2_biomass_pixel <- phase1_bleach_biomass_pixel
names(phase2_biomass_pixel) <- "phase2_biomass_g_per_pixel"  # grams per pixel

# Scaling factor for Phase 2 (mid-term restoration goal)
S2 <- 6.9

# 1) Masks and definitions
# baseline_degraded: pixels that were degraded in the baseline (class == 1)
baseline_degraded <- (baseline_class == 1)

# improved_phase2: baseline-degraded pixels that are healthy in Phase 2
# This covers:
#  - pixels that recovered during Phase 2 (baseline degraded -> Phase2 healthy)
#  - pixels restored during Phase 1 that remained healthy into Phase 2
improved_phase2 <- baseline_degraded & (phase2_class == 2)

# 2) Ensure always-healthy pixels retain baseline biomass
# Define pixels that are healthy across all scenario rasters (baseline, phase1, bleach, phase2)
always_healthy <- (baseline_class == 2) & (phase1_class == 2) & (bleach_class == 2) & (phase2_class == 2)

# Overwrite those always-healthy pixels with the baseline per-pixel biomass
# (baseline_biomass_pixel should have been created earlier as g per pixel)
# Robust enforcement of "always healthy keep baseline" rule

if (!exists("baseline_biomass_pixel")) {
  warning("baseline_biomass_pixel not found: cannot enforce 'always healthy keep baseline' rule.")
} else {
  # Ensure baseline_biomass_pixel is same geometry as phase2 (if not, resample)
  same_geom <- terra::compareGeom(baseline_biomass_pixel, phase2_biomass_pixel, stopOnError = FALSE)
  if (!same_geom) {
    message("baseline_biomass_pixel and phase2_biomass_pixel differ in geometry. Resampling baseline to phase2 grid (nearest).")
    baseline_for_use <- terra::resample(baseline_biomass_pixel, phase2_biomass_pixel, method = "near")
  } else {
    baseline_for_use <- baseline_biomass_pixel
  }
  
  # Only overwrite where baseline has non-NA values and pixel is always_healthy
  overwrite_mask <- always_healthy & (!is.na(baseline_for_use))
  
  # Count pixels to be overwritten (quick diagnostic)
  n_overwrite <- tryCatch({
    as.integer(global(overwrite_mask, "sum", na.rm = TRUE)[1])
  }, error = function(e) NA_integer_)
  
  message("Always-healthy pixels (total): ", tryCatch(as.integer(global(always_healthy, "sum", na.rm = TRUE)[1]), error = function(e) NA_integer_))
  message("Pixels to be overwritten with baseline biomass: ", n_overwrite)
  
  # Do the overwrite
  phase2_biomass_pixel[overwrite_mask] <- baseline_for_use[overwrite_mask]
  
  # Optional: warn if many always_healthy pixels still NA in baseline
  n_always_healthy <- tryCatch(as.integer(global(always_healthy, "sum", na.rm = TRUE)[1]), error = function(e) NA_integer_)
  if (!is.na(n_always_healthy) && !is.na(n_overwrite) && n_overwrite < n_always_healthy) {
    message("Note: some always-healthy pixels lacked baseline per-pixel biomass (left as NA). Consider checking baseline_biomass_pixel coverage.")
  }
}

# 3) Apply S2 scaling (6.9×) to improved pixels by zone
# For each zone, use the mean degraded biomass (g/100 m^2) -> convert to g/pixel, then * S2
for (z in all_zones) {
  if (is.na(z)) next
  
  zone_row <- zone_table %>% filter(Zone == z)
  if (nrow(zone_row) == 0) next
  zid <- zone_row$zone_id
  
  zone_pixels <- zones_raster == zid
  improved_zone <- zone_pixels & improved_phase2
  
  # mean degraded biomass for this zone (g / 100 m^2)
  bm_deg <- biomass_lookup %>% filter(Zone == z, class == 1) %>% pull(biomass)
  if (length(bm_deg) == 0) bm_deg <- NA
  
  # convert to grams per pixel (using pixel area conversion function defined earlier)
  bm_deg_g_per_pixel <- if (!is.na(bm_deg)) g100m2_to_g_per_pixel(bm_deg) else NA_real_
  
  # Assign S2 * degraded baseline biomass to improved pixels
  phase2_biomass_pixel[improved_zone] <- bm_deg_g_per_pixel * S2
}

# Depth rule:
# deep > 10 m: +12% biomass
# 4) Depth-based adjustment for pixels degraded by bleaching and not restored by Phase 2
# deep_mask already defined in section 8c above

# Pixels that were healthy in Phase 1 but degraded by bleaching (phase1 -> bleach)
bleach_degraded <- (phase1_class == 2) & (bleach_class == 1)

# Of those, pixels that remained degraded after Phase 2
bleach_degraded_not_restored <- bleach_degraded & (phase2_class == 1)

# Deep recovery pixels: bleached-degraded, not restored, and deep
deep_recovery_pixels <- bleach_degraded_not_restored & deep_mask

# Apply +12% to deep recovery pixels relative to post-bleach biomass
phase2_biomass_pixel[deep_recovery_pixels] <-
  phase1_bleach_biomass_pixel[deep_recovery_pixels] * 1.12

# Note: shallow pixels that are bleach_degraded_not_restored keep their post-bleach values
# because phase2_biomass_pixel was initially copied from phase1_bleach_biomass_pixel.

# 5) Quick diagnostics (recommended)
n_cells <- ncell(phase2_biomass_pixel)
na_cells <- sum(is.na(values(phase2_biomass_pixel)))
message("Phase2: total cells = ", n_cells, "; NA cells = ", na_cells)

total_g_phase2 <- tryCatch(global(phase2_biomass_pixel, fun = "sum", na.rm = TRUE)[1], error = function(e) NA_real_)
if (!is.na(total_g_phase2)) {
  total_tons_phase2 <- total_g_phase2 / 1e6     # g -> metric tons
  message("Phase2 total biomass: ", signif(total_tons_phase2, 6), " metric tons")
} else {
  message("Could not compute global sum for phase2_biomass_pixel.")
}

# Save phase2 raster
writeRaster(phase2_biomass_pixel,
            filename = file.path(spatial_out_dir, "phase2_fish_biomass_g_per_pixel.tif"),
            overwrite = TRUE)

# ELIMINAR
phase2_biomass_pixel <- rast("C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/spatial_files/phase2_fish_biomass_g100m2.tif")

############################################################
# 10. BLOCK A: TOTAL biomass per zone (tons)
############################################################

# call the safe_zonal function to sum all values from baseline_biomass_pixel within each fishing zone defined by zones_raster
baseline_total <- safe_zonal(baseline_biomass_pixel, 
                             zones_raster, 
                             "baseline_g_sum", 
                             fun_type = "sum") %>%
  mutate(baseline_kg = g_to_kg(baseline_g_sum), 
         baseline_ton = kg_to_ton(baseline_kg)
         )

# similar sum but for phase 1 restoration scenario
phase1_total <- safe_zonal(phase1_biomass_pixel, zones_raster, "phase1_g_sum", fun_type = "sum") %>%
  mutate(phase1_kg = g_to_kg(phase1_g_sum), phase1_ton = kg_to_ton(phase1_kg))

# sum after bleaching event
phase1_bleach_total <- safe_zonal(phase1_bleach_biomass_pixel, 
                                  zones_raster, 
                                  "phase1_bleach_g_sum", 
                                  fun_type = "sum") %>%
  mutate(phase1_bleach_kg = g_to_kg(phase1_bleach_g_sum), 
         phase1_bleach_ton = kg_to_ton(phase1_bleach_kg)
         )

# sum for bleaching_only scenario (no restoration, bleaching on baseline)
bleaching_only_total <- safe_zonal(bleaching_only_biomass_pixel, 
                                   zones_raster, 
                                   "bleaching_only_g_sum", 
                                   fun_type = "sum") %>%
  mutate(bleaching_only_kg = g_to_kg(bleaching_only_g_sum), 
         bleaching_only_ton = kg_to_ton(bleaching_only_kg)
         )

# sum for bleaching_only ~5 Years scenario (deep +12% partial recovery, no coral restoration)
bleaching_only_phase2_total <- safe_zonal(bleaching_only_phase2_biomass_pixel, 
                                          zones_raster, 
                                          "bleaching_only_phase2_g_sum", 
                                          fun_type = "sum") %>%
  mutate(bleaching_only_phase2_kg = g_to_kg(bleaching_only_phase2_g_sum), 
         bleaching_only_phase2_ton = kg_to_ton(bleaching_only_phase2_kg)
         )

# sum after phase 2 restoration
phase2_total <- safe_zonal(phase2_biomass_pixel, 
                           zones_raster, 
                           "phase2_g_sum", 
                           fun_type = "sum") %>%
  mutate(phase2_kg = g_to_kg(phase2_g_sum), 
         phase2_ton = kg_to_ton(phase2_kg)
         )

# create a table will al the results
biomass_total_compare <- baseline_total %>%
  left_join(phase1_total %>% dplyr::select(zone_id, phase1_ton), by = "zone_id") %>%
  left_join(bleaching_only_total %>% dplyr::select(zone_id, bleaching_only_ton), by = "zone_id") %>%
  left_join(bleaching_only_phase2_total %>% dplyr::select(zone_id, bleaching_only_phase2_ton), by = "zone_id") %>%
  left_join(phase1_bleach_total %>% dplyr::select(zone_id, phase1_bleach_ton), by = "zone_id") %>%
  left_join(phase2_total %>% dplyr::select(zone_id, phase2_ton), by = "zone_id") %>%
  left_join(zone_table, by = "zone_id") %>%
  dplyr::select(Zone, zone_id, baseline_ton, phase1_ton, bleaching_only_ton, bleaching_only_phase2_ton, phase1_bleach_ton, phase2_ton)

write_xlsx(biomass_total_compare, path = file.path(tabplot_out_dir, "fish_biomass_total_by_zone_tons.xlsx"))

# Plot totals
total_long <- biomass_total_compare %>%
  pivot_longer(
    cols = c(baseline_ton, phase1_ton, bleaching_only_ton, bleaching_only_phase2_ton, phase1_bleach_ton, phase2_ton),
    names_to = "Scenario",
    values_to = "Biomass_ton"
  ) %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = c("baseline_ton", "phase1_ton", "bleaching_only_ton", "bleaching_only_phase2_ton", "phase1_bleach_ton", "phase2_ton"),
      labels = c("Baseline", "Phase 1", "Bleaching Only", "Bleaching Only Phase 2", "Bleaching", "Phase 2")
    )
  )

write.csv(total_long, file.path(out_dir, "total_biomass_zone_scenarios.csv"), row.names = FALSE)

## Two-panel comparative bar plot: Bleaching Only vs With Restoration
total_long_2panel <- bind_rows(
  total_long %>%
    filter(Scenario %in% c("Baseline", "Bleaching Only")) %>%
    mutate(scenario_group = "No Intervention"),
  total_long %>%
    filter(Scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(scenario_group = "With Restoration")
) %>%
  mutate(scenario_group = factor(scenario_group, levels = c("No Intervention", "With Restoration")))

### Two-panel line plot: Bleaching Only vs With Restoration
# Temporal x-axis labels mirror those used in the coral-coverage analysis:
#   Baseline = Year 0 | 1 Year | 2 Years | ~5 Years
#
# Bleaching Only panel — 4 synthetic time points:
#   "Baseline"  = baseline values (Year 0)
#   "1 Year"    = baseline values (no change without restoration)
#   "2 Years"   = bleaching_only values (bleaching impact on baseline)
#   "~5 Years"  = bleaching_only_phase2 values (deep +12% recovery, no coral restoration)
#
# With Restoration panel — same values, renamed to temporal labels:
#   "Baseline"  = Baseline | "1 Year" = Phase 1 | "2 Years" = Bleaching | "~5 Years" = Phase 2

library(scales)
library(dplyr)
library(ggplot2)
library(forcats)

# Zone colour palette (Dark2) — derived from total_long zone levels
zone_levels <- levels(factor(total_long$Zone))
if (is.null(zone_levels)) zone_levels <- sort(unique(as.character(total_long$Zone)))
zone_levels <- as.character(zone_levels)
zone_cols <- setNames(
  RColorBrewer::brewer.pal(n = length(zone_levels), name = "Dark2"),
  zone_levels
)

# Build two-panel long data with temporal x-axis labels
total_long2_2panel <- bind_rows(

  # ---- Bleaching Only panel (4 synthetic time points) ----
  total_long %>%
    filter(Scenario == "Baseline") %>%
    mutate(Scenario = "Year zero", scenario_group = "No Intervention"),
  total_long %>%
    filter(Scenario == "Baseline") %>%
    mutate(Scenario = "1 Year", scenario_group = "No Intervention"),
  total_long %>%
    filter(Scenario == "Bleaching Only") %>%
    mutate(Scenario = "2 Years", scenario_group = "No Intervention"),
  total_long %>%
    filter(Scenario == "Bleaching Only Phase 2") %>%
    mutate(Scenario = "~5 Years", scenario_group = "No Intervention"),

  # ---- With Restoration panel (renamed to temporal labels) ----
  total_long %>%
    filter(Scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(
      Scenario = recode(as.character(Scenario),
                        "Baseline" = "Year zero",
                        "Phase 1"  = "1 Year",
                        "Bleaching" = "2 Years",
                        "Phase 2"  = "~5 Years"),
      scenario_group = "With Restoration"
    )

) %>%
  mutate(
    Scenario       = factor(Scenario, levels = c("Year zero", "1 Year", "2 Years", "~5 Years")),
    scenario_group = factor(scenario_group, levels = c("No Intervention", "With Restoration")),
    Zone           = factor(Zone)
  )

### National Biomass Sum per scenario group (dashed overlay line)
national_biomass_2panel <- total_long2_2panel %>%
  group_by(scenario_group, Scenario) %>%
  summarise(Biomass_ton = sum(Biomass_ton, na.rm = TRUE), .groups = "drop") %>%
  mutate(Zone = "National")

p_total_line_2panel <- ggplot(total_long2_2panel,
                              # Cambiamos a factor directamente
                              aes(x = Scenario, y = Biomass_ton, group = Zone, color = Zone)) +
  
  # Capas regionales
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +
  
  # Capa Nacional (Dashed black line)
  # Usamos group = 1 para que ggplot entienda que debe conectar los puntos del factor
  geom_line(data = national_biomass_2panel,
            aes(x = Scenario, y = Biomass_ton, group = 1, linetype = "National total"),
            linewidth = 1.2, color = "black") +
  geom_point(data = national_biomass_2panel,
             aes(x = Scenario, y = Biomass_ton),
             size = 3, shape = 18, color = "black") +
  
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  
  # Escala Y (Asegúrate de tener cargada la librería 'scales' para label_number)
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  
  scale_color_manual(values = zone_cols) +
  scale_linetype_manual(
    name   = NULL,
    values = c("National total" = "dashed"),
    guide  = guide_legend(
      override.aes = list(color = "black", linewidth = 1.2, shape = NA)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Biomass (tons)",
    color = "Zone"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 13)
  )

p_total_line_2panel

# Output folder for publication figures
out_dir_fig <- "G:/Shared drives/NSF CoPE internal/2 - Deliverables/Publications/Olaya_et_al_Belize_FisheryModel/figures"

ggsave(
  filename = file.path(out_dir, "Fig_biomass_tons_byFishing_zone.svg"),
  plot = p_total_line_2panel,
  device = svglite::svglite,
  width = 12,
  height = 5.2,
  units = "in",
  fix_text_size = FALSE
)

