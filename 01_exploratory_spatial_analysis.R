# ============================================================
# Script: 01_exploratory_spatial_analysis.R
# Purpose: Exploratory spatial analysis of terrain and coral
#          reef datasets for St. Thomas (USVI).
#
# Study area : St. Thomas, US Virgin Islands (USVI)
# Target CRS : WGS 84 / UTM Zone 20N  (EPSG:32620)
# Resolution : 50 m × 50 m (working resolution)
#
# Tasks
# ------
# 1. Clip the global WCMC coral-reef mask to the extent of the
#    bathymetry raster and save the result as `coral_mask_stt`.
# 2. Report the resolution and coordinate reference system (CRS)
#    of every raster (.tif) and vector (.shp) input file.
# 3. Count how many hard-coral survey points fall inside the
#    clipped coral mask and compute the proportion.
#
# Author : [your name]
# Date   : 2024
# ============================================================


# ── 0. Required packages ──────────────────────────────────────────────────────
# terra  – modern raster / vector handling in R
# sf     – simple-features vector analysis
# dplyr  – tidy data manipulation (optional convenience helpers)

if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("sf",    quietly = TRUE)) install.packages("sf")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(terra)   # SpatRaster / SpatVector
library(sf)      # Simple Features
library(dplyr)   # Data wrangling


# ── 1. File paths ─────────────────────────────────────────────────────────────
# Adjust the root path to match the location of your shared drive.
# All paths below follow the folder structure described in the project brief.

root <- "G:/Shared drives/NSF CoPE internal/GIS_CoPE/GIS_USVI"

# ---- Raster inputs -----------------------------------------------------------
path_bathy          <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "00_bathymetry_source/STTSTJ_2m.tif")

path_slope_10m      <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "00_intermediate/slope_10m.tif")

path_slope_10m_buf  <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "00_intermediate/slope_10m_240m_buffer.tif")

path_slope_30m      <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "00_intermediate/slope_30m.tif")

path_slope_50m      <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "01_slope/slope_50m.tif")

path_slope_240m_50m <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "01_slope/slope_240m_50m.tif")

path_curv_2m        <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "02_curvature/curvature_2m.tif")

path_slope_of_slope <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "02_curvature/slope_of_slope_50m.tif")

path_aspect_2m      <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "03_aspect/aspect_2m.tif")

# ---- Vector inputs -----------------------------------------------------------
path_coral_wcmc <- file.path(root,
  "0_source_data_usvi/WCMC008_CoralReefs2021_v4_1",
  "14_001_WCMC008_CoralReefs2021_v4_1/01_Data",
  "WCMC008_CoralReef2021_Py_v4_1.shp")

path_hard_corals <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/StThomas_Benthic",
  "hard_corals_data_filtered.shp")


# ── 2. Load data ──────────────────────────────────────────────────────────────
message("Loading raster and vector data …")

# Load the high-resolution (2 m) bathymetry raster.
# This raster defines the spatial footprint of the study area.
bathy <- rast(path_bathy)

# Load all remaining terrain rasters.
slope_10m      <- rast(path_slope_10m)
slope_10m_buf  <- rast(path_slope_10m_buf)
slope_30m      <- rast(path_slope_30m)
slope_50m      <- rast(path_slope_50m)
slope_240m_50m <- rast(path_slope_240m_50m)
curv_2m        <- rast(path_curv_2m)
slope_of_slope <- rast(path_slope_of_slope)
aspect_2m      <- rast(path_aspect_2m)

# Load vector layers.
# `st_read` imports shapefiles as sf data frames.
coral_wcmc   <- st_read(path_coral_wcmc,   quiet = TRUE)
hard_corals  <- st_read(path_hard_corals,  quiet = TRUE)

message("All files loaded successfully.")


# ── TASK 1 ──────────────────────────────────────────────────────────────────
# Clip the global WCMC coral-reef polygon mask to the bounding box
# of the bathymetry raster (i.e., the study area extent).
# Result: coral_mask_stt  – coral polygons restricted to St. Thomas / STJ.
# ─────────────────────────────────────────────────────────────────────────────

message("\n── TASK 1: Clipping coral mask to bathymetry extent ──────────────────")

# Step 1a: Reproject the WCMC coral polygons to the target CRS (EPSG:32620)
#          so that all spatial operations are performed in the same system.
target_crs <- "EPSG:32620"   # WGS 84 / UTM Zone 20N

coral_wcmc_proj <- st_transform(coral_wcmc, crs = target_crs)

# Step 1b: Extract the bounding box of the bathymetry raster and convert it
#          to an sf polygon so it can be used as a clipping geometry.
bathy_bbox    <- st_bbox(project(bathy, target_crs))  # reproject raster extent
bathy_extent  <- st_as_sfc(bathy_bbox)                # bbox → sfc polygon
st_crs(bathy_extent) <- target_crs

# Step 1c: Clip (intersect) the coral polygons with the bathymetry extent.
#          Only coral polygons (or polygon parts) that fall within the
#          bathymetry footprint are retained.
coral_mask_stt <- st_intersection(coral_wcmc_proj, bathy_extent)

message(paste("coral_mask_stt: ", nrow(coral_mask_stt), "polygon(s) after clipping."))

# Step 1d: Optionally save the clipped mask to disk.
#          Uncomment the line below to write the result as a shapefile.
# st_write(coral_mask_stt, "coral_mask_stt.shp", delete_dsn = TRUE)


# ── TASK 2 ──────────────────────────────────────────────────────────────────
# Verify the spatial resolution and CRS of all raster and vector files.
# This is an essential quality-control step before any analysis to ensure
# that all layers are consistent (or to identify layers that need reprojection
# / resampling before use together).
# ─────────────────────────────────────────────────────────────────────────────

message("\n── TASK 2: CRS and resolution check ──────────────────────────────────")

# Helper function: extract key metadata from a terra SpatRaster.
raster_info <- function(r, label) {
  res_xy <- res(r)
  data.frame(
    file       = label,
    type       = "raster",
    crs_name   = crs(r, describe = TRUE)$name,
    epsg       = crs(r, describe = TRUE)$code,
    res_x_m    = round(res_xy[1], 4),
    res_y_m    = round(res_xy[2], 4),
    nrow       = nrow(r),
    ncol       = ncol(r),
    stringsAsFactors = FALSE
  )
}

# Helper function: extract key metadata from an sf object.
vector_info <- function(v, label) {
  crs_obj  <- st_crs(v)
  data.frame(
    file       = label,
    type       = "vector",
    crs_name   = crs_obj$Name,
    epsg       = crs_obj$epsg,
    res_x_m    = NA_real_,
    res_y_m    = NA_real_,
    nrow       = nrow(v),
    ncol       = NA_integer_,
    stringsAsFactors = FALSE
  )
}

# Build a summary table for all rasters.
raster_summary <- dplyr::bind_rows(
  raster_info(bathy,          "STTSTJ_2m.tif"),
  raster_info(slope_10m,      "slope_10m.tif"),
  raster_info(slope_10m_buf,  "slope_10m_240m_buffer.tif"),
  raster_info(slope_30m,      "slope_30m.tif"),
  raster_info(slope_50m,      "slope_50m.tif"),
  raster_info(slope_240m_50m, "slope_240m_50m.tif"),
  raster_info(curv_2m,        "curvature_2m.tif"),
  raster_info(slope_of_slope, "slope_of_slope_50m.tif"),
  raster_info(aspect_2m,      "aspect_2m.tif")
)

# Build a summary table for all vector layers
# (including the freshly clipped coral_mask_stt).
vector_summary <- dplyr::bind_rows(
  vector_info(coral_wcmc,      "WCMC008_CoralReef2021_Py_v4_1.shp (original)"),
  vector_info(coral_mask_stt,  "coral_mask_stt (clipped to STT)"),
  vector_info(hard_corals,     "hard_corals_data_filtered.shp")
)

# Display the summaries.
message("\nRaster summary:")
print(raster_summary)

message("\nVector summary:")
print(vector_summary)

# Flag any rasters that are NOT already in the target CRS.
non_target_rasters <- raster_summary[
  !is.na(raster_summary$epsg) & raster_summary$epsg != "32620", ]

if (nrow(non_target_rasters) > 0) {
  message("\nWARNING – The following rasters are NOT in EPSG:32620:")
  print(non_target_rasters[, c("file", "crs_name", "epsg")])
} else {
  message("\nAll rasters are in the target CRS (EPSG:32620).")
}

# Flag rasters that do NOT have the 50 m × 50 m working resolution.
non_50m_rasters <- raster_summary[
  raster_summary$res_x_m != 50 | raster_summary$res_y_m != 50, ]

if (nrow(non_50m_rasters) > 0) {
  message("\nNOTE – The following rasters differ from the 50 m working resolution:")
  print(non_50m_rasters[, c("file", "res_x_m", "res_y_m")])
}


# ── TASK 3 ──────────────────────────────────────────────────────────────────
# Count hard-coral survey points that fall INSIDE the clipped coral mask
# (coral_mask_stt) and compute the proportion relative to all points.
#
# Rationale: The WCMC coral-reef polygon layer represents the known spatial
# distribution of coral reefs.  Survey points that fall within these polygons
# are considered to be on confirmed reef habitat.  Points outside may represent
# adjacent or cryptic reef habitats, data entry errors, or off-reef surveys.
# ─────────────────────────────────────────────────────────────────────────────

message("\n── TASK 3: Survey points vs. coral mask ──────────────────────────────")

# Step 3a: Ensure the survey points are in the same CRS as coral_mask_stt.
hard_corals_proj <- st_transform(hard_corals, crs = target_crs)

# Step 3b: Perform a spatial join.
#          st_join with join = st_within returns only the rows from
#          hard_corals_proj for which the geometry lies within a polygon
#          of coral_mask_stt.  We use left = TRUE so that ALL points are
#          kept; those outside the mask will have NA in the joined columns.
points_joined <- st_join(
  hard_corals_proj,
  coral_mask_stt,
  join      = st_within,
  left      = TRUE,     # keep ALL survey points
  suffix    = c("_pt", "_mask")
)

# Step 3c: A survey point is considered "inside the coral mask" when the
#          spatial join produced at least one matching polygon row
#          (i.e., the joined attribute from coral_mask_stt is NOT NA).
#          We use any geometry column from coral_mask_stt as indicator.
#          If coral_mask_stt has no unique ID column, fall back to checking
#          duplicate row indices introduced by the join.

# Identify a column that comes from coral_mask_stt (not the point layer).
# After st_join the new columns have suffix "_mask"; pick the first one.
mask_cols <- grep("_mask$", names(points_joined), value = TRUE)

if (length(mask_cols) > 0) {
  # Use the first mask-origin column: NA means the point was outside.
  inside_flag <- !is.na(points_joined[[mask_cols[1]]])
} else {
  # Fallback: duplicate rows after st_join indicate overlap with a polygon.
  # Points that matched ≥1 polygon appear more than once; use unique row index.
  inside_flag <- !duplicated(st_drop_geometry(hard_corals_proj))
  warning("Could not detect mask columns; falling back to row-duplicate method.")
}

n_total   <- nrow(hard_corals_proj)
n_inside  <- sum(inside_flag, na.rm = TRUE)
n_outside <- n_total - n_inside
prop_inside <- n_inside / n_total

# Step 3d: Report the results.
message(sprintf("\nTotal hard-coral survey points  : %d", n_total))
message(sprintf("Points INSIDE  coral mask (stt) : %d", n_inside))
message(sprintf("Points OUTSIDE coral mask (stt) : %d", n_outside))
message(sprintf("Proportion inside coral mask    : %.1f%%  (%.4f)",
                prop_inside * 100, prop_inside))

# Step 3e: Create a summary data frame for further use or export.
results_summary <- data.frame(
  category        = c("Inside coral mask", "Outside coral mask", "Total"),
  n_points        = c(n_inside, n_outside, n_total),
  proportion      = c(prop_inside, 1 - prop_inside, 1)
)

message("\nSummary table:")
print(results_summary)


# ── End of script ─────────────────────────────────────────────────────────────
message("\nExploratory spatial analysis complete.")
