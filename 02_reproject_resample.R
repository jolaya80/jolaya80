# ============================================================
# Script: 02_reproject_resample.R
# Purpose: Reproject all terrain rasters from NAD83 / UTM Zone 20N
#          (EPSG:26920) to WGS 84 / UTM Zone 20N (EPSG:32620), and
#          resample them from 2 m native resolution to the 50 m × 50 m
#          working resolution required for the NOAA benthic survey model.
#
# Study area : St. Thomas, US Virgin Islands (USVI)
# Source CRS : NAD83 / UTM Zone 20N  (EPSG:26920)  ← current state
# Target CRS : WGS 84 / UTM Zone 20N (EPSG:32620)  ← required
# Native res : 2 m × 2 m                           ← current state
# Target res : 50 m × 50 m                         ← required
#
# ── WHY ARE THE TWO CRS DIFFERENT? ───────────────────────────────────────────
# NAD83 (EPSG:26920) and WGS 84 (EPSG:32620) share the same UTM Zone 20N
# projection and are nearly identical for practical purposes in the Caribbean
# (datum shift < 1 m).  However, they carry different datum definitions, so
# any analysis that mixes both will generate a CRS-mismatch warning and may
# produce subtle alignment errors at the sub-metre level.  Reprojecting to a
# single CRS (WGS 84 / UTM 20N) is the correct approach.
#
# ── RISKS OF REPROJECTION ─────────────────────────────────────────────────────
# 1. Datum shift (very small here, < 1 m):
#    NAD83 ↔ WGS 84 offsets in the USVI are negligible relative to the 50 m
#    working resolution, so no meaningful spatial error is introduced.
#
# 2. Pixel re-alignment during reprojection:
#    Even when the datum shift is sub-metre, `project()` must interpolate new
#    cell values because the output grid corners shift slightly.  Using
#    "bilinear" interpolation for continuous variables (depth, slope, curvature,
#    aspect) is standard and accurate; "near" (nearest-neighbour) should be
#    used for categorical rasters (none here).
#
# 3. Resolution change (2 m → 50 m): *** MOST IMPORTANT RISK ***
#    Aggregating a 2 m raster to 50 m means each output cell summarises a
#    25 × 25 = 625 source cells.  The aggregation method matters:
#
#    Variable          Recommended method   Reason
#    ──────────────    ─────────────────    ──────────────────────────────────
#    Bathymetry        mean                 Average depth of the 50 m cell
#    Slope (any)       mean                 Average slope within the cell
#    Curvature         mean                 Average curvature
#    Slope-of-slope    mean                 Average second-order derivative
#    Aspect            circular mean*       Aspect is a circular (angular)
#                                           variable; a standard mean would
#                                           be incorrect near 0°/360°.
#                                           (*) terra does not have a built-in
#                                           circular mean aggregation; a
#                                           sin/cos decomposition workaround
#                                           is implemented below.
#
#    Using "mean" for aspect would give nonsensical results where aspect
#    crosses the 0°/360° boundary (e.g., averaging 1° and 359° gives 180°
#    instead of 0°).  The circular mean fix below converts to unit vectors,
#    aggregates, and converts back.
#
# 4. Loss of fine-scale variability:
#    Smoothing from 2 m to 50 m irreversibly removes high-frequency bottom
#    topography detail.  This is intentional and necessary for the 50 m model,
#    but the original high-resolution files should always be kept as archive.
#    RECOMMENDATION: write outputs to a *new* folder; never overwrite originals.
#
# 5. Increased file size:
#    Reprojection + resampling actually REDUCES file size (625× fewer cells),
#    but the intermediate reprojected-at-2m grid can be large in memory.
#    The workflow below reprojects and aggregates in a single pass where
#    possible to minimise peak RAM usage.
#
# ── WORKFLOW ──────────────────────────────────────────────────────────────────
# For each raster:
#   Step A – project() : reproject from EPSG:26920 → EPSG:32620 (keeps 2 m res)
#   Step B – aggregate(): resample from 2 m → 50 m using the appropriate method
#   Step C – writeRaster(): save the result to the output folder
#
# ============================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
library(terra)


# ── 1. Paths ──────────────────────────────────────────────────────────────────
root <- "G:/Shared drives/NSF CoPE internal/GIS_CoPE/GIS_USVI"

# ---- Input rasters (all at 2 m, EPSG:26920) ----------------------------------
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

# ---- Output folder -----------------------------------------------------------
# Outputs are saved here; originals are NEVER overwritten.
out_dir <- file.path(root,
  "2_model_inputs_usvi/Fisheries/coral_cover_modeling",
  "05_preparation_spatial_predictors/02_terrain_analysis_outputs",
  "50m_WGS84_outputs")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message("Output folder: ", out_dir)


# ── 2. Constants ──────────────────────────────────────────────────────────────
target_crs <- "EPSG:32620"   # WGS 84 / UTM Zone 20N
native_res  <- 2             # original pixel size in metres
target_res  <- 50            # desired pixel size in metres
agg_factor  <- target_res / native_res   # = 25 (each output cell = 25×25 input cells)

message(sprintf(
  "Aggregation factor: %d (each 50 m cell summarises %d × %d = %d source cells)",
  agg_factor, agg_factor, agg_factor, agg_factor^2))


# ── 3. Helper functions ────────────────────────────────────────────────────────

# reproject_resample()
# --------------------
# Reprojects a SpatRaster from its native CRS (EPSG:26920) to the target CRS
# (EPSG:32620) and then aggregates it to 50 m using the specified method.
#
# Arguments:
#   r       – input SpatRaster
#   method  – aggregation method: "mean", "median", "min", "max", or "modal"
#   label   – short name used in progress messages
#
# Returns a SpatRaster at 50 m / EPSG:32620.
reproject_resample <- function(r, method = "mean", label = "raster") {

  message(sprintf("  [%s] Reprojecting from %s → %s …",
                  label,
                  crs(r, describe = TRUE)$code,
                  target_crs))

  # Step A: Reproject at native resolution (2 m).
  # Using "bilinear" interpolation for continuous surfaces.
  r_proj <- project(r, target_crs, method = "bilinear")

  message(sprintf("  [%s] Aggregating %d m → %d m using '%s' …",
                  label, native_res, target_res, method))

  # Step B: Aggregate from 2 m to 50 m.
  # fact = 25 means 25 cells in each direction are combined per output cell.
  # na.rm = TRUE: output cell is computed even when some source cells are NA
  #               (e.g., at the raster boundary).
  r_50m <- aggregate(r_proj, fact = agg_factor, fun = method, na.rm = TRUE)

  message(sprintf("  [%s] Done. Output resolution: %s × %s m, CRS: %s",
                  label,
                  res(r_50m)[1], res(r_50m)[2],
                  crs(r_50m, describe = TRUE)$code))

  return(r_50m)
}


# circular_mean_aspect()
# ----------------------
# Aspect is a *circular* (angular) variable in degrees [0, 360).
# A naive arithmetic mean fails near the 0°/360° wrap-around:
#   mean(1°, 359°) = 180°  ← wrong; correct answer is 0° (North).
#
# The circular mean is computed by:
#   1. Converting each degree value to unit-vector components (sin, cos).
#   2. Averaging the sin and cos components independently (standard mean).
#   3. Recovering the angle with atan2(mean_sin, mean_cos).
#
# This function performs that operation via terra's built-in raster maths.
circular_mean_aspect <- function(r_aspect, label = "aspect") {

  message(sprintf("  [%s] Reprojecting aspect …", label))
  r_proj <- project(r_aspect, target_crs, method = "bilinear")

  # Convert degrees to radians for trigonometric functions.
  r_rad <- r_proj * (pi / 180)

  # Decompose into unit-vector components.
  r_sin <- sin(r_rad)
  r_cos <- cos(r_rad)

  message(sprintf("  [%s] Aggregating aspect via circular mean …", label))

  # Aggregate the sin and cos planes independently using arithmetic mean.
  sin_agg <- aggregate(r_sin, fact = agg_factor, fun = "mean", na.rm = TRUE)
  cos_agg <- aggregate(r_cos, fact = agg_factor, fun = "mean", na.rm = TRUE)

  # Recover the mean angle in degrees [0, 360).
  # atan2 returns values in (−π, π]; adding 2π and taking modulo 2π
  # maps the result to [0, 2π), then convert back to degrees.
  r_mean_rad <- (atan2(sin_agg, cos_agg) + 2 * pi) %% (2 * pi)
  r_mean_deg <- r_mean_rad * (180 / pi)

  message(sprintf("  [%s] Done. Circular-mean aspect at 50 m computed.", label))

  return(r_mean_deg)
}


# ── 4. Process each raster ────────────────────────────────────────────────────
# All terrain variables use mean aggregation EXCEPT aspect, which requires
# the circular mean to handle the 0°/360° wrap-around correctly.

message("\n── Processing STTSTJ_2m (bathymetry) ─────────────────────────────────")
bathy_50m <- reproject_resample(rast(path_bathy), method = "mean", label = "bathy")
writeRaster(bathy_50m, file.path(out_dir, "bathymetry_50m.tif"), overwrite = TRUE)

message("\n── Processing slope_10m ──────────────────────────────────────────────")
slope_10m_50m <- reproject_resample(rast(path_slope_10m), method = "mean", label = "slope_10m")
writeRaster(slope_10m_50m, file.path(out_dir, "slope_10m_50m.tif"), overwrite = TRUE)

message("\n── Processing slope_10m_240m_buffer ──────────────────────────────────")
slope_10m_buf_50m <- reproject_resample(rast(path_slope_10m_buf), method = "mean",
                                        label = "slope_10m_buf")
writeRaster(slope_10m_buf_50m, file.path(out_dir, "slope_10m_240m_buffer_50m.tif"),
            overwrite = TRUE)

message("\n── Processing slope_30m ──────────────────────────────────────────────")
slope_30m_50m <- reproject_resample(rast(path_slope_30m), method = "mean", label = "slope_30m")
writeRaster(slope_30m_50m, file.path(out_dir, "slope_30m_50m.tif"), overwrite = TRUE)

message("\n── Processing slope_50m ──────────────────────────────────────────────")
# NOTE: slope_50m was *computed* at a 50 m analysis window but stored at 2 m
# native resolution; we still aggregate to get true 50 m output cells.
slope_50m_50m <- reproject_resample(rast(path_slope_50m), method = "mean", label = "slope_50m")
writeRaster(slope_50m_50m, file.path(out_dir, "slope_50m_50m.tif"), overwrite = TRUE)

message("\n── Processing slope_240m_50m ─────────────────────────────────────────")
slope_240m_50m_50m <- reproject_resample(rast(path_slope_240m_50m), method = "mean",
                                         label = "slope_240m")
writeRaster(slope_240m_50m_50m, file.path(out_dir, "slope_240m_50m_50m.tif"),
            overwrite = TRUE)

message("\n── Processing curvature_2m ───────────────────────────────────────────")
# Curvature can be negative (concave) or positive (convex); mean is appropriate.
curv_50m <- reproject_resample(rast(path_curv_2m), method = "mean", label = "curvature")
writeRaster(curv_50m, file.path(out_dir, "curvature_50m.tif"), overwrite = TRUE)

message("\n── Processing slope_of_slope_50m ─────────────────────────────────────")
sos_50m <- reproject_resample(rast(path_slope_of_slope), method = "mean",
                               label = "slope_of_slope")
writeRaster(sos_50m, file.path(out_dir, "slope_of_slope_50m_50m.tif"), overwrite = TRUE)

message("\n── Processing aspect_2m (circular mean) ──────────────────────────────")
# Aspect uses circular-mean aggregation; see circular_mean_aspect() above.
aspect_50m <- circular_mean_aspect(rast(path_aspect_2m), label = "aspect")
writeRaster(aspect_50m, file.path(out_dir, "aspect_50m.tif"), overwrite = TRUE)


# ── 5. Verification ───────────────────────────────────────────────────────────
# After processing, verify that all output rasters share the correct CRS and
# resolution.  Any deviation indicates an error in the pipeline.

message("\n── Verification: checking output rasters ─────────────────────────────")

output_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)

verify_raster <- function(path) {
  r <- rast(path)
  data.frame(
    file     = basename(path),
    epsg     = crs(r, describe = TRUE)$code,
    res_x_m  = round(res(r)[1], 2),
    res_y_m  = round(res(r)[2], 2),
    nrow     = nrow(r),
    ncol     = ncol(r),
    stringsAsFactors = FALSE
  )
}

verification_table <- do.call(rbind, lapply(output_files, verify_raster))
print(verification_table)

# Pass / fail check.
all_ok <- all(
  verification_table$epsg    == "32620" &
  verification_table$res_x_m == target_res &
  verification_table$res_y_m == target_res
)

if (all_ok) {
  message("\n✓ All output rasters are in EPSG:32620 at 50 m × 50 m resolution.")
} else {
  warning("\n✗ One or more output rasters do NOT meet the target CRS/resolution.")
  print(verification_table[
    verification_table$epsg != "32620" |
    verification_table$res_x_m != target_res |
    verification_table$res_y_m != target_res, ])
}

message("\nReprojection and resampling complete.")
message("Output folder: ", out_dir)
