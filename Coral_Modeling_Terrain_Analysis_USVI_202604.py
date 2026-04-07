"""
Coral Modeling – Terrain Analysis, USVI
========================================
St. Thomas / St. John bathymetric terrain predictor stack

VERSION NOTES (April 2026 – corrected)
---------------------------------------
All fixes recommended in the April 2026 code review are applied:

Fix 1  – WorkflowConfig: single-step aggregation factor (2m→50m = ×25);
          correct 240m focal radius at 2m resolution (120 cells, not 5);
          remove confusing CELL_FACTOR_2M_TO_6M / CELL_FACTOR_6M_TO_50M.

Fix 2  – BLOCK 1: reproject bathymetry to target CRS (EPSG:32620) before
          any terrain analysis.  All derivatives inherit the correct datum.

Fix 3  – BLOCK 3: single Aggregate(factor=25) step produces a *true* 50m
          slope raster.  The previous two-step chain (×3 copy, then a
          BlockStatistics moving-window that does NOT change cell size)
          left the output at ~6m or 2m resolution with a wrong filename.

Fix 4  – BLOCK 4: NbrCircle radius changed from 5 cells to 120 cells so
          the focal neighbourhood covers the intended 240m on a 2m raster.

Fix 5  – BLOCK 5: slope_of_slope computed from the corrected slope_50m.tif
          (truly 50m pixels), not from a 2m copy labelled "50m".

Fix 6  – BLOCK 6 & 6b: aspect computed at 2m (correct), then aggregated to
          50m using circular-mean decomposition (sin/cos) so the 0°/360°
          wrap-around is handled correctly.

Fix 7  – BLOCK 7 & 7b: curvature computed at 2m (correct), then aggregated
          to 50m with standard mean.

Output naming convention
-------------------------
  *_2m.tif   – computed at native 2m resolution (intermediate or archive)
  *_50m.tif  – final product; every cell is a true 50m × 50m pixel in
               EPSG:32620 (WGS 84 / UTM Zone 20N)

Dependencies
------------
  ArcGIS Pro  + Spatial Analyst Extension
  arcpy (bundled with ArcGIS Pro)
  pandas (optional, for tabular summaries)

BLOCK EXECUTION ORDER
----------------------
  BLOCK 0  → BLOCK 1 → BLOCK 2 → VERIFICATION 2A
  → BLOCK 3 → BLOCK 4 → BLOCK 5 → BLOCK 6 → BLOCK 6b
  → BLOCK 7 → BLOCK 7b → BLOCK 8
"""

import arcpy
import os
import sys
import math
from datetime import datetime


# ============================================================================
# LOGGING UTILITY
# ============================================================================

class Logger:
    """Timestamped console logging."""

    @staticmethod
    def info(message):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] ✓ {message}")

    @staticmethod
    def warning(message):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] ⚠ {message}")

    @staticmethod
    def error(message):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] ✗ {message}")

    @staticmethod
    def header(message):
        print("\n" + "=" * 70)
        print(f"  {message}")
        print("=" * 70)


# ============================================================================
# BLOCK 0 – CONFIGURATION  (Fix 1)
# ============================================================================
"""
BLOCK 0: CONFIGURATION
======================
Run this FIRST — before all other blocks.

FIX 1 applied
-------------
• WORKING_RESOLUTION target is 50m; source is 2m → aggregation factor = 25.
• FOCAL_RADIUS_240M is expressed in cells AT THE 2m RASTER, not at 50m.
  On the 2m raster: 240m / 2m = 120 cells.
• CELL_FACTOR_2M_TO_6M and CELL_FACTOR_6M_TO_50M removed (caused the
  wrong two-step aggregation that left outputs at 6m or 2m).
"""

class WorkflowConfig:
    """Configuration for St. Thomas benthic habitat analysis."""

    # ── Study Area ──────────────────────────────────────────────────────────
    SITE_NAME          = "St_Thomas_USVI"

    # FIX 1 & 2: target CRS for ALL outputs
    TARGET_CRS_EPSG    = 32620          # WGS 84 / UTM Zone 20N
    TARGET_CRS_STR     = "EPSG:32620"

    # FIX 1: resolutions
    SOURCE_RESOLUTION  = 2    # native bathymetry pixel size (metres)
    TARGET_RESOLUTION  = 50   # final output pixel size (metres)
    AGG_FACTOR         = 25   # 50 / 2 = 25  (each output cell = 25×25 source cells)

    # FIX 1: focal radii expressed in cells ON THE 2m SOURCE RASTER
    # Neighbourhood of 30m  → 30  / 2 = 15 cells
    # Neighbourhood of 240m → 240 / 2 = 120 cells
    FOCAL_RADIUS_30M_at2m  = 15   # cells on 2m raster
    FOCAL_RADIUS_240M_at2m = 120  # cells on 2m raster  ← was incorrectly 5

    # ── Project Directories ─────────────────────────────────────────────────
    PROJECT_ROOT = (
        r"G:\Shared drives\NSF CoPE internal\GIS_CoPE\GIS_USVI"
        r"\2_model_inputs_usvi\Fisheries\coral_cover_modeling"
        r"\05_preparation_spatial_predictors"
    )
    OUTPUTS_ROOT = os.path.join(PROJECT_ROOT, "02_terrain_analysis_outputs")

    SOURCE_DATA_ROOT = (
        r"G:\Shared drives\NSF CoPE internal\GIS_CoPE\GIS_USVI\0_source_data_usvi"
    )

    # ── Bathymetry Source (direct TIF — no GDB export needed) ───────────────
    BATHYMETRY_SOURCE_TIF = (
        r"G:\Shared drives\NSF CoPE internal\GIS_CoPE\GIS_USVI"
        r"\2_model_inputs_usvi\Fisheries\coral_cover_modeling"
        r"\05_preparation_spatial_predictors\02_terrain_analysis_outputs"
        r"\00_bathymetry_source\STTSTJ_2m.tif"
    )

    # BATHYMETRY_SOURCE_2M is kept for downstream compatibility; it points
    # directly to the source TIF (no intermediate copy required).
    BATHYMETRY_WORKING_DIR = os.path.join(OUTPUTS_ROOT, "00_bathymetry_source")
    BATHYMETRY_SOURCE_2M   = BATHYMETRY_SOURCE_TIF

    # FIX 2: reprojected bathymetry in EPSG:32620 — used as snap raster and
    #         input for ALL terrain derivatives
    BATHYMETRY_WGS84_2M    = os.path.join(BATHYMETRY_WORKING_DIR, "STTSTJ_2m_EPSG32620.tif")

    # ── Optional Data Sources ───────────────────────────────────────────────
    TSS_SOURCE      = None
    WAVE_ERA5_POINTS = None
    AOI_BOUNDARY    = None

    # ── IDW parameters (wave interpolation) ─────────────────────────────────
    IDW_POWER            = 2
    IDW_MAX_NEIGHBORS    = 4
    IDW_SEMIMAJOR_AXIS   = 150000   # metres
    IDW_SEMIMINOR_AXIS   = 75000    # metres


# Instantiate global config
config = WorkflowConfig()

arcpy.env.overwriteOutput = True

try:
    arcpy.CheckExtension("Spatial")
    Logger.info("Spatial Analyst Extension available")
except Exception:
    Logger.warning("Spatial Analyst Extension may not be available")

Logger.header("BLOCK 0: CONFIGURATION LOADED")
Logger.info(f"Site             : {config.SITE_NAME}")
Logger.info(f"Target CRS       : {config.TARGET_CRS_STR}")
Logger.info(f"Source resolution: {config.SOURCE_RESOLUTION}m")
Logger.info(f"Target resolution: {config.TARGET_RESOLUTION}m")
Logger.info(f"Aggregation factor: {config.AGG_FACTOR}")
Logger.info(f"240m focal radius : {config.FOCAL_RADIUS_240M_at2m} cells "
            f"({config.FOCAL_RADIUS_240M_at2m * config.SOURCE_RESOLUTION}m)")
Logger.info("\n✓ Configuration ready. Proceed to BLOCK 1.")


# ============================================================================
# BLOCK 1 – INITIALIZE WORKSPACE & EXPORT / REPROJECT BATHYMETRY  (Fix 2)
# ============================================================================
"""
BLOCK 1: INITIALIZE WORKSPACE & EXPORT / REPROJECT BATHYMETRY
==============================================================
Steps
-----
1. Create output directory tree.
2. Export the STTSTJ_2m raster from the GDB to a working .tif
   (native CRS, 2m — kept as archive).
3. FIX 2: Reproject the exported .tif to EPSG:32620 (WGS 84 / UTM 20N)
   if the native CRS is NAD83 / UTM 20N (EPSG:26920).
   NAD83 ↔ WGS84 datum shift is < 1m in the USVI, but mixing the two
   CRS labels causes alignment warnings in downstream tools.
4. Set snap raster to the reprojected 2m bathymetry.

All subsequent terrain calculations use BATHYMETRY_WGS84_2M.
"""

Logger.header("BLOCK 1: INITIALIZE WORKSPACE & EXPORT / REPROJECT BATHYMETRY")

# ── 1. Create output directories ─────────────────────────────────────────────
Logger.info("Creating output directory structure...")

output_dirs = {
    'bathymetry'   : config.BATHYMETRY_WORKING_DIR,
    'intermediate' : os.path.join(config.OUTPUTS_ROOT, "00_intermediate"),
    'slope'        : os.path.join(config.OUTPUTS_ROOT, "01_slope"),
    'curvature'    : os.path.join(config.OUTPUTS_ROOT, "02_curvature"),
    'aspect'       : os.path.join(config.OUTPUTS_ROOT, "03_aspect"),
    'wave'         : os.path.join(config.OUTPUTS_ROOT, "04_wave"),
    'tss'          : os.path.join(config.OUTPUTS_ROOT, "05_tss"),
}

for key, dir_path in output_dirs.items():
    os.makedirs(dir_path, exist_ok=True)
    Logger.info(f"  ✓ {key}: {dir_path}")

# ── 2. Verify bathymetry source TIF exists ────────────────────────────────────
if os.path.exists(config.BATHYMETRY_SOURCE_2M):
    Logger.info(f"✓ Source bathymetry .tif found: {config.BATHYMETRY_SOURCE_2M}")
else:
    Logger.error(f"Source bathymetry not found: {config.BATHYMETRY_SOURCE_2M}")
    raise FileNotFoundError(
        f"STTSTJ_2m.tif not found at expected path:\n  {config.BATHYMETRY_SOURCE_2M}"
    )

# ── 3. FIX 2 – Reproject to EPSG:32620 if needed ─────────────────────────────
if os.path.exists(config.BATHYMETRY_WGS84_2M):
    Logger.info("✓ Reprojected bathymetry .tif already exists (cached)")
else:
    Logger.info("Checking CRS of exported bathymetry…")
    bathy_desc = arcpy.Describe(config.BATHYMETRY_SOURCE_2M)
    native_epsg = (bathy_desc.spatialReference.factoryCode
                   if bathy_desc.spatialReference else None)
    Logger.info(f"  Native EPSG: {native_epsg}")

    if native_epsg == config.TARGET_CRS_EPSG:
        Logger.info("  CRS already matches target — creating symlink-equivalent copy")
        arcpy.management.CopyRaster(config.BATHYMETRY_SOURCE_2M,
                                    config.BATHYMETRY_WGS84_2M)
    else:
        Logger.warning(
            f"  CRS mismatch (native={native_epsg}, target={config.TARGET_CRS_EPSG}) "
            f"— reprojecting…"
        )
        # Project raster: bilinear for continuous bathymetry surface
        arcpy.management.ProjectRaster(
            in_raster           = config.BATHYMETRY_SOURCE_2M,
            out_raster          = config.BATHYMETRY_WGS84_2M,
            out_coor_system     = arcpy.SpatialReference(config.TARGET_CRS_EPSG),
            resampling_type     = "BILINEAR",
            cell_size           = str(config.SOURCE_RESOLUTION),  # keep 2m
        )
        Logger.info(f"✓ Reprojected to EPSG:{config.TARGET_CRS_EPSG}")

# ── 4. Set workspace and snap raster ─────────────────────────────────────────
arcpy.env.workspace    = config.OUTPUTS_ROOT
arcpy.env.snapRaster   = config.BATHYMETRY_WGS84_2M
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(config.TARGET_CRS_EPSG)

Logger.info(f"✓ Snap raster : {os.path.basename(config.BATHYMETRY_WGS84_2M)}")
Logger.info(f"✓ Output CRS  : EPSG:{config.TARGET_CRS_EPSG}")

# ── 5. Verify ─────────────────────────────────────────────────────────────────
try:
    desc    = arcpy.Describe(config.BATHYMETRY_WGS84_2M)
    raster  = arcpy.Raster(config.BATHYMETRY_WGS84_2M)
    Logger.info(f"  Cell size  : {raster.meanCellWidth:.1f}m × {raster.meanCellHeight:.1f}m")
    Logger.info(f"  CRS        : {desc.spatialReference.name}")
    Logger.info(f"  Depth range: {float(raster.minimum):.2f}m to {float(raster.maximum):.2f}m")
except Exception as e:
    Logger.warning(f"Could not read full raster properties: {e}")

Logger.header("✓ BLOCK 1 COMPLETE")
Logger.info("Bathymetry ready at 2m / EPSG:32620.")
Logger.info("Proceed to BLOCK 2.")


# ============================================================================
# BLOCK 2 – VALIDATE INPUTS
# ============================================================================
"""
BLOCK 2: VALIDATE INPUTS & CHECK DATA SOURCES
==============================================
Confirms which processing steps can proceed.
"""

Logger.header("BLOCK 2: VALIDATE INPUTS & CHECK DATA SOURCES")

data_status = {'bathymetry': 'ERROR', 'tss': 'NOT_FOUND',
               'wave': 'NOT_FOUND', 'aoi': 'NOT_FOUND'}

if os.path.exists(config.BATHYMETRY_WGS84_2M):
    Logger.info("✓ BATHYMETRY (EPSG:32620, 2m): Ready")
    data_status['bathymetry'] = 'READY'
else:
    Logger.error(f"✗ BATHYMETRY not found: {config.BATHYMETRY_WGS84_2M}")
    raise FileNotFoundError("Bathymetry required to proceed")

for key, src, label in [
    ('tss',  config.TSS_SOURCE,       "TSS"),
    ('wave', config.WAVE_ERA5_POINTS, "WAVE ERA5 points"),
    ('aoi',  config.AOI_BOUNDARY,     "AOI boundary"),
]:
    if src and os.path.exists(src):
        Logger.info(f"✓ {label}: Ready")
        data_status[key] = 'READY'
    else:
        Logger.warning(f"⚠ {label}: Not found (step will be skipped)")
        data_status[key] = 'SKIPPED'

Logger.header("✓ BLOCK 2 COMPLETE")
Logger.info("Proceed to BLOCK 3.")


# ============================================================================
# VERIFICATION 2A – BATHYMETRY INSPECTION
# ============================================================================

Logger.header("VERIFICATION 2A: BATHYMETRY INSPECTION")

try:
    desc   = arcpy.Describe(config.BATHYMETRY_WGS84_2M)
    raster = arcpy.Raster(config.BATHYMETRY_WGS84_2M)
    print(f"\n  File      : {os.path.basename(config.BATHYMETRY_WGS84_2M)}")
    print(f"  CRS       : {desc.spatialReference.name}")
    print(f"  EPSG      : {desc.spatialReference.factoryCode}")
    print(f"  Cell size : {raster.meanCellWidth:.1f}m × {raster.meanCellHeight:.1f}m")
    print(f"  Depth min : {float(raster.minimum):.2f}m")
    print(f"  Depth max : {float(raster.maximum):.2f}m")
    print(f"  Depth mean: {float(raster.mean):.2f}m")
    print(f"  Depth std : {float(raster.standardDeviation):.2f}m")
    print(f"  Extent    : {desc.extent}")
except Exception as e:
    Logger.warning(f"Minor issue reading properties: {e}")

Logger.header("✓ VERIFICATION 2A COMPLETE")
Logger.info("Proceed to BLOCK 3.")


# ============================================================================
# HELPER: convert negative bathymetry to positive once (used by Slope/Aspect/Curvature)
# ============================================================================

_bathy_positive_path = os.path.join(output_dirs['intermediate'],
                                    "STTSTJ_2m_positive_EPSG32620.tif")

def get_positive_bathymetry():
    """Return path to positive-value bathymetry (depth→elevation flip), creating if needed."""
    if not os.path.exists(_bathy_positive_path):
        Logger.info("Converting bathymetry to positive values (depth→elevation)…")
        from arcpy.sa import Raster as SAR
        bathy = SAR(config.BATHYMETRY_WGS84_2M)
        positive = bathy * -1
        positive.save(_bathy_positive_path)
        Logger.info(f"✓ Positive bathymetry saved: {_bathy_positive_path}")
    return _bathy_positive_path


# ============================================================================
# BLOCK 3 – SLOPE AT MULTIPLE SCALES  (Fix 3)
# ============================================================================
"""
BLOCK 3: CALCULATE SLOPE AT MULTIPLE SCALES
============================================
FIX 3 applied
-------------
The previous code produced:
  slope_10m.tif  → actually 2m  (Slope_3d on 2m DEM → 2m output)
  slope_30m.tif  → actually 2m  (BlockStatistics moving-window; does NOT
                                  change cell size)
  slope_50m.tif  → actually 2m  (Copy of slope_30m)

The corrected approach:
  1. Compute slope at native 2m resolution   → slope_2m.tif   (intermediate)
  2. Aggregate 2m → 50m in ONE step (×25)    → slope_50m.tif  (final, true 50m)

Output files
  00_intermediate/slope_2m.tif       – slope at native 2m resolution
  01_slope/slope_50m.tif             – mean slope aggregated to true 50m

PROCESSING TIME: ~60-90 seconds
"""

from arcpy.sa import (Raster, FocalStatistics, NbrCircle, NbrRectangle,
                      BlockStatistics, Slope, Aspect, Curvature)

Logger.header("BLOCK 3: CALCULATE SLOPE AT MULTIPLE SCALES")

bathymetry_positive = get_positive_bathymetry()

# ── Step 1: slope at native 2m resolution ─────────────────────────────────────
slope_2m = os.path.join(output_dirs['intermediate'], "slope_2m.tif")

Logger.info("[1/2] Computing slope at 2m resolution…")
Logger.info("  ⏳ ~30-40 seconds…")
try:
    arcpy.env.snapRaster = config.BATHYMETRY_WGS84_2M
    arcpy.Slope_3d(bathymetry_positive, slope_2m,
                   output_measurement="DEGREE", z_factor=1.0)
    Logger.info(f"✓ slope_2m.tif created ({os.path.getsize(slope_2m)/1e6:.1f} MB)")
    slope_rast = arcpy.Raster(slope_2m)
    Logger.info(f"  Slope range: {float(slope_rast.minimum):.2f}° – "
                f"{float(slope_rast.maximum):.2f}°")
except Exception as e:
    Logger.error(f"Slope_3d failed: {e}")
    raise

# ── Step 2: aggregate 2m → 50m (factor = 25, mean) ───────────────────────────
# FIX 3: single Aggregate step; output cells are truly 50m × 50m.
slope_50m = os.path.join(output_dirs['slope'], "slope_50m.tif")

Logger.info(f"\n[2/2] Aggregating slope 2m → 50m (factor={config.AGG_FACTOR}, mean)…")
Logger.info("  ⏳ ~20-30 seconds…")
try:
    from arcpy.sa import Aggregate as SA_Aggregate
    result = SA_Aggregate(Raster(slope_2m), config.AGG_FACTOR, "MEAN",
                          extent_handling="EXPAND", ignore_nodata="DATA")
    result.save(slope_50m)
    Logger.info(f"✓ slope_50m.tif created ({os.path.getsize(slope_50m)/1e6:.1f} MB)")
    # Verify true 50m cell size
    d = arcpy.Describe(slope_50m)
    _r = arcpy.Raster(slope_50m)
    Logger.info(f"  Cell size : {_r.meanCellWidth:.1f}m × {_r.meanCellHeight:.1f}m")
    Logger.info(f"  CRS       : {d.spatialReference.name}")
except Exception as e:
    Logger.error(f"Aggregation failed: {e}")
    raise

Logger.header("✓ BLOCK 3 COMPLETE")
Logger.info("Outputs:")
Logger.info(f"  00_intermediate/slope_2m.tif   (2m, EPSG:32620)")
Logger.info(f"  01_slope/slope_50m.tif         (50m, EPSG:32620, true 50m cells)")
Logger.info("Next: BLOCK 4 (Slope at 240m scale)")


# ============================================================================
# BLOCK 4 – SLOPE AT 240M SCALE  (Fix 4)
# ============================================================================
"""
BLOCK 4: CALCULATE SLOPE AT 240M SCALE
=======================================
FIX 4 applied
-------------
The previous code used NbrCircle(5, "CELL") on the 2m raster:
  5 cells × 2m = 10m radius → 20m neighbourhood  (NOT 240m)

Correct radius for a 240m neighbourhood on a 2m raster:
  240m / 2m = 120 cells  → FOCAL_RADIUS_240M_at2m = 120

This block:
  1. Applies FocalStatistics with NbrCircle(120) on slope_2m.tif.
  2. Aggregates the focal-smoothed 2m raster to 50m (factor=25, mean).

Output files
  00_intermediate/slope_2m_240m_focal.tif  – focal mean at 240m scale, 2m pixels
  01_slope/slope_240m_50m.tif              – aggregated to true 50m output

PROCESSING TIME: ~5-8 min (large focal radius on a 2m raster)
"""

Logger.header("BLOCK 4: CALCULATE SLOPE AT 240M SCALE")

slope_2m_240m_focal = os.path.join(output_dirs['intermediate'],
                                   "slope_2m_240m_focal.tif")
slope_240m_50m      = os.path.join(output_dirs['slope'], "slope_240m_50m.tif")

# ── Step 1: focal mean with 240m radius on 2m raster ─────────────────────────
Logger.info(f"[1/2] FocalStatistics (radius={config.FOCAL_RADIUS_240M_at2m} cells "
            f"= {config.FOCAL_RADIUS_240M_at2m * config.SOURCE_RESOLUTION}m) "
            f"on slope_2m.tif…")
Logger.warning("  ⏳ Large focal radius — expect 5-8 minutes, please wait…")
try:
    arcpy.env.snapRaster = config.BATHYMETRY_WGS84_2M
    neighborhood = NbrCircle(config.FOCAL_RADIUS_240M_at2m, "CELL")  # FIX 4
    focal_result = FocalStatistics(Raster(slope_2m), neighborhood, "MEAN")
    focal_result.save(slope_2m_240m_focal)
    Logger.info(f"✓ 240m focal result saved ({os.path.getsize(slope_2m_240m_focal)/1e6:.1f} MB)")
except Exception as e:
    Logger.error(f"FocalStatistics failed: {e}")
    raise

# ── Step 2: aggregate 2m → 50m ────────────────────────────────────────────────
Logger.info(f"\n[2/2] Aggregating 2m → 50m (factor={config.AGG_FACTOR})…")
try:
    result = SA_Aggregate(Raster(slope_2m_240m_focal), config.AGG_FACTOR, "MEAN",
                          extent_handling="EXPAND", ignore_nodata="DATA")
    result.save(slope_240m_50m)
    Logger.info(f"✓ slope_240m_50m.tif created ({os.path.getsize(slope_240m_50m)/1e6:.1f} MB)")
    d = arcpy.Describe(slope_240m_50m)
    Logger.info(f"  Cell size : {arcpy.Raster(slope_240m_50m).meanCellWidth:.1f}m × {arcpy.Raster(slope_240m_50m).meanCellHeight:.1f}m")
except Exception as e:
    Logger.error(f"Aggregation failed: {e}")
    raise

Logger.header("✓ BLOCK 4 COMPLETE")
Logger.info("Outputs:")
Logger.info(f"  00_intermediate/slope_2m_240m_focal.tif  (2m, EPSG:32620)")
Logger.info(f"  01_slope/slope_240m_50m.tif              (50m, EPSG:32620, true 50m)")
Logger.info("Next: BLOCK 5 (Terrain Ruggedness)")


# ============================================================================
# BLOCK 5 – SLOPE OF SLOPE (TERRAIN RUGGEDNESS)  (Fix 5)
# ============================================================================
"""
BLOCK 5: SLOPE OF SLOPE (TERRAIN RUGGEDNESS)
=============================================
FIX 5 applied
-------------
The previous code computed Slope() on slope_50m.tif when that file was
actually a 2m copy labelled "50m", yielding near-zero slope-of-slope values
(gradient of a nearly-flat smoothed surface).

The corrected approach:
  • Use the TRUE 50m slope (slope_50m.tif from FIX 3).
  • Slope() on a true 50m raster captures meaningful rate-of-change in slope.
  • Output is already 50m — no further aggregation needed.

Output files
  02_curvature/slope_of_slope_50m.tif  (50m, EPSG:32620)

PROCESSING TIME: ~30-40 seconds
"""

Logger.header("BLOCK 5: SLOPE OF SLOPE (TERRAIN RUGGEDNESS)")

slope_of_slope_50m = os.path.join(output_dirs['curvature'], "slope_of_slope_50m.tif")

Logger.info("[1/1] Computing slope of slope from slope_50m.tif (true 50m)…")
Logger.info("  ⏳ ~30-40 seconds…")
try:
    arcpy.env.snapRaster = config.BATHYMETRY_WGS84_2M
    # FIX 5: input is the TRUE 50m slope produced in BLOCK 3
    ruggedness = Slope(Raster(slope_50m), z_factor=1.0)
    ruggedness.save(slope_of_slope_50m)
    Logger.info(f"✓ slope_of_slope_50m.tif created "
                f"({os.path.getsize(slope_of_slope_50m)/1e6:.1f} MB)")
    d = arcpy.Describe(slope_of_slope_50m)
    Logger.info(f"  Cell size : {arcpy.Raster(slope_of_slope_50m).meanCellWidth:.1f}m × {arcpy.Raster(slope_of_slope_50m).meanCellHeight:.1f}m")
except Exception as e:
    Logger.error(f"Failed: {e}")
    raise

Logger.header("✓ BLOCK 5 COMPLETE")
Logger.info("Output:")
Logger.info(f"  02_curvature/slope_of_slope_50m.tif  (50m, EPSG:32620)")
Logger.info("Next: BLOCK 6 (Aspect)")


# ============================================================================
# BLOCK 6 – ASPECT AT 2M  +  BLOCK 6b – AGGREGATE TO 50M (circular mean)
#                                                           (Fix 6)
# ============================================================================
"""
BLOCK 6: ASPECT (DIRECTIONAL EXPOSURE) — 2m native resolution
BLOCK 6b: AGGREGATE ASPECT TO 50m — circular mean via sin/cos decomposition

FIX 6 applied
-------------
Aspect is a CIRCULAR variable (0–360°, wraps at 0°/360°).
A naive arithmetic mean is wrong near the wrap-around:
  mean(1°, 359°) = 180°  ← incorrect; correct answer is 0° (North).

Correct approach (circular mean):
  1. Decompose aspect into unit-vector components:
       aspect_sin = sin(aspect_rad)
       aspect_cos = cos(aspect_rad)
  2. Aggregate sin and cos components independently (standard mean).
  3. Recover the circular mean angle:
       mean_aspect = atan2(mean_sin, mean_cos) converted to [0°, 360°).

Output files
  03_aspect/aspect_2m.tif    (2m, EPSG:32620)        ← archive
  03_aspect/aspect_50m.tif   (50m, EPSG:32620, circular mean)  ← final model input
"""

import math

Logger.header("BLOCK 6 + 6b: ASPECT — 2m THEN AGGREGATED TO 50m (circular mean)")

aspect_2m  = os.path.join(output_dirs['aspect'], "aspect_2m.tif")
aspect_sin_2m = os.path.join(output_dirs['intermediate'], "aspect_sin_2m.tif")
aspect_cos_2m = os.path.join(output_dirs['intermediate'], "aspect_cos_2m.tif")
aspect_50m = os.path.join(output_dirs['aspect'], "aspect_50m.tif")

# ── 6: Compute aspect at 2m ───────────────────────────────────────────────────
Logger.info("[1/4] Computing aspect at 2m resolution…")
Logger.info("  ⏳ ~30-40 seconds…")
try:
    arcpy.env.snapRaster = config.BATHYMETRY_WGS84_2M
    aspect_result = Aspect(Raster(bathymetry_positive))
    aspect_result.save(aspect_2m)
    Logger.info(f"✓ aspect_2m.tif created ({os.path.getsize(aspect_2m)/1e6:.1f} MB)")
    d = arcpy.Describe(aspect_2m)
    Logger.info(f"  Cell size : {arcpy.Raster(aspect_2m).meanCellWidth:.1f}m × {arcpy.Raster(aspect_2m).meanCellHeight:.1f}m")
except Exception as e:
    Logger.error(f"Aspect failed: {e}")
    raise

# ── 6b: Circular mean aggregation ─────────────────────────────────────────────
# Step 1: Convert aspect (degrees) → radians → sin / cos components
Logger.info("\n[2/4] Decomposing aspect into sin/cos components…")
try:
    asp = Raster(aspect_2m)

    # The Aspect tool returns -1 for flat areas; mask those out
    from arcpy.sa import SetNull, IsNull
    asp_valid = SetNull(asp < 0, asp)               # mask flat-area pixels

    DEG2RAD = math.pi / 180.0
    asp_rad = asp_valid * DEG2RAD
    sin_r = arcpy.sa.Sin(asp_rad)
    cos_r = arcpy.sa.Cos(asp_rad)
    sin_r.save(aspect_sin_2m)
    cos_r.save(aspect_cos_2m)
    Logger.info("✓ sin/cos components saved")
except Exception as e:
    Logger.error(f"Sin/cos decomposition failed: {e}")
    raise

# Step 2: Aggregate sin and cos independently (standard mean)
Logger.info("\n[3/4] Aggregating sin/cos components 2m → 50m…")
sin_50m_path = os.path.join(output_dirs['intermediate'], "aspect_sin_50m.tif")
cos_50m_path = os.path.join(output_dirs['intermediate'], "aspect_cos_50m.tif")
try:
    sin_50 = SA_Aggregate(Raster(aspect_sin_2m), config.AGG_FACTOR, "MEAN",
                          extent_handling="EXPAND", ignore_nodata="DATA")
    sin_50.save(sin_50m_path)
    cos_50 = SA_Aggregate(Raster(aspect_cos_2m), config.AGG_FACTOR, "MEAN",
                          extent_handling="EXPAND", ignore_nodata="DATA")
    cos_50.save(cos_50m_path)
    Logger.info("✓ sin/cos aggregated to 50m")
except Exception as e:
    Logger.error(f"Sin/cos aggregation failed: {e}")
    raise

# Step 3: Recover circular mean angle: atan2(sin, cos) → [0°, 360°)
Logger.info("\n[4/4] Recovering circular-mean aspect angle [0°, 360°)…")
try:
    sin50 = Raster(sin_50m_path)
    cos50 = Raster(cos_50m_path)

    RAD2DEG = 180.0 / math.pi
    # atan2 returns (−π, π]; shift to [0, 2π) then convert to degrees
    mean_rad = arcpy.sa.ATan2(sin50, cos50)              # radians, −π to π
    TWO_PI   = 2.0 * math.pi
    mean_deg = ((mean_rad + TWO_PI) % TWO_PI) * RAD2DEG  # degrees [0, 360)
    mean_deg.save(aspect_50m)
    Logger.info(f"✓ aspect_50m.tif created ({os.path.getsize(aspect_50m)/1e6:.1f} MB)")
    d = arcpy.Describe(aspect_50m)
    _r = arcpy.Raster(aspect_50m)
    Logger.info(f"  Cell size : {_r.meanCellWidth:.1f}m × {_r.meanCellHeight:.1f}m")
    Logger.info(f"  CRS       : {d.spatialReference.name}")
except Exception as e:
    Logger.error(f"Circular mean recovery failed: {e}")
    raise

Logger.header("✓ BLOCK 6 + 6b COMPLETE")
Logger.info("Outputs:")
Logger.info(f"  03_aspect/aspect_2m.tif   (2m, EPSG:32620, archive)")
Logger.info(f"  03_aspect/aspect_50m.tif  (50m, EPSG:32620, circular mean)")
Logger.info("Next: BLOCK 7 (Curvature)")


# ============================================================================
# BLOCK 7 – CURVATURE AT 2M  +  BLOCK 7b – AGGREGATE TO 50M  (Fix 7)
# ============================================================================
"""
BLOCK 7: CURVATURE — 2m native resolution
BLOCK 7b: AGGREGATE CURVATURE TO 50m (standard mean)

FIX 7 applied
-------------
Previously only curvature_2m.tif was produced; the R reproject script was
relied upon to aggregate it to 50m.  Now aggregation is done here in Python
so every final model input is 50m / EPSG:32620 without needing the R script.

ECOLOGICAL INTERPRETATION
  Positive curvature = convex (ridges, peaks)  → flow / sediment divergence
  Negative curvature = concave (valleys)       → flow / sediment convergence
  Zero curvature     = flat / linear terrain

Output files
  02_curvature/curvature_2m.tif   (2m, EPSG:32620)   ← archive
  02_curvature/curvature_50m.tif  (50m, EPSG:32620)  ← final model input

PROCESSING TIME: ~30-40 seconds (2m), ~20 seconds (aggregate)
"""

Logger.header("BLOCK 7 + 7b: CURVATURE — 2m THEN AGGREGATED TO 50m")

curvature_2m  = os.path.join(output_dirs['curvature'], "curvature_2m.tif")
curvature_50m = os.path.join(output_dirs['curvature'], "curvature_50m.tif")

# ── 7: Compute curvature at 2m ────────────────────────────────────────────────
Logger.info("[1/2] Computing curvature at 2m resolution…")
Logger.info("  ⏳ ~30-40 seconds…")
try:
    arcpy.env.snapRaster = config.BATHYMETRY_WGS84_2M
    curv_result = Curvature(Raster(bathymetry_positive))
    curv_result.save(curvature_2m)
    Logger.info(f"✓ curvature_2m.tif created ({os.path.getsize(curvature_2m)/1e6:.1f} MB)")
    d = arcpy.Describe(curvature_2m)
    Logger.info(f"  Cell size : {arcpy.Raster(curvature_2m).meanCellWidth:.1f}m × {arcpy.Raster(curvature_2m).meanCellHeight:.1f}m")
except Exception as e:
    Logger.error(f"Curvature failed: {e}")
    raise

# ── 7b: Aggregate 2m → 50m (standard mean; curvature can be negative) ─────────
Logger.info("\n[2/2] Aggregating curvature 2m → 50m (factor=25, mean)…")
try:
    result = SA_Aggregate(Raster(curvature_2m), config.AGG_FACTOR, "MEAN",
                          extent_handling="EXPAND", ignore_nodata="DATA")
    result.save(curvature_50m)
    Logger.info(f"✓ curvature_50m.tif created ({os.path.getsize(curvature_50m)/1e6:.1f} MB)")
    d = arcpy.Describe(curvature_50m)
    _r = arcpy.Raster(curvature_50m)
    Logger.info(f"  Cell size : {_r.meanCellWidth:.1f}m × {_r.meanCellHeight:.1f}m")
    Logger.info(f"  CRS       : {d.spatialReference.name}")
except Exception as e:
    Logger.error(f"Aggregation failed: {e}")
    raise

Logger.header("✓ BLOCK 7 + 7b COMPLETE")
Logger.info("Outputs:")
Logger.info(f"  02_curvature/curvature_2m.tif   (2m, EPSG:32620, archive)")
Logger.info(f"  02_curvature/curvature_50m.tif  (50m, EPSG:32620)")
Logger.info("Next: BLOCK 8 (Final verification)")


# ============================================================================
# BLOCK 8 – FINAL VERIFICATION
# ============================================================================
"""
BLOCK 8: FINAL VERIFICATION
============================
Reports CRS and cell size for every final 50m output.
All files should show:
  EPSG: 32620
  Cell size: 50.0m × 50.0m
"""

Logger.header("BLOCK 8: FINAL VERIFICATION — ALL 50m OUTPUTS")

final_outputs = {
    "slope_50m.tif"          : slope_50m,
    "slope_240m_50m.tif"     : slope_240m_50m,
    "slope_of_slope_50m.tif" : slope_of_slope_50m,
    "aspect_50m.tif"         : aspect_50m,
    "curvature_50m.tif"      : curvature_50m,
}

print(f"\n{'File':<30} {'CRS / EPSG':<35} {'Cell X (m)':>10} {'Cell Y (m)':>10}")
print("-" * 90)

all_ok = True
for label, path in final_outputs.items():
    if not os.path.exists(path):
        print(f"  {'MISSING':<28} {label}")
        all_ok = False
        continue
    try:
        d    = arcpy.Describe(path)
        _r   = arcpy.Raster(path)
        epsg = d.spatialReference.factoryCode if d.spatialReference else "?"
        cx   = _r.meanCellWidth
        cy   = _r.meanCellHeight
        crs_label = f"{d.spatialReference.name} ({epsg})"
        ok = (epsg == config.TARGET_CRS_EPSG
              and abs(cx - config.TARGET_RESOLUTION) < 1
              and abs(cy - config.TARGET_RESOLUTION) < 1)
        flag = "✓" if ok else "✗"
        print(f"  {flag} {label:<28} {crs_label:<35} {cx:>10.1f} {cy:>10.1f}")
        if not ok:
            all_ok = False
    except Exception as e:
        print(f"  ? {label:<28} Error: {e}")
        all_ok = False

print()
if all_ok:
    Logger.info("ALL outputs are in EPSG:32620 at 50m × 50m.")
    Logger.info("The R script 02_reproject_resample.R is no longer needed.")
else:
    Logger.warning("One or more outputs have CRS or resolution issues — check above.")

Logger.header("WORKFLOW COMPLETE")
Logger.info("Final 50m / EPSG:32620 model inputs are ready.")
Logger.info("")
Logger.info("Intermediate 2m archives (for reference / QA):")
Logger.info(f"  00_bathymetry_source/STTSTJ_2m_native.tif")
Logger.info(f"  00_bathymetry_source/STTSTJ_2m_EPSG32620.tif")
Logger.info(f"  00_intermediate/slope_2m.tif")
Logger.info(f"  00_intermediate/slope_2m_240m_focal.tif")
Logger.info(f"  03_aspect/aspect_2m.tif")
Logger.info(f"  02_curvature/curvature_2m.tif")
