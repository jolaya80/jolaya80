# Reproducibility Report
## Coral Reef Terrain Analysis — St. Thomas / St. John, US Virgin Islands (USVI)
### Spatial Predictor Stack Construction at 50 m / EPSG:32620

---

**Authors:** J. Olaya (domain scientist) · GitHub Copilot Coding Agent (AI assistant)  
**Date:** April 2026  
**Repository:** [jolaya80/jolaya80](https://github.com/jolaya80/jolaya80)  
**Primary script:** `Coral_Modeling_Terrain_Analysis_USVI_202604.py`  
**Supporting scripts:** `01_exploratory_spatial_analysis.R`, `02_reproject_resample.R`

---

## Table of Contents
1. [Purpose of This Document](#1-purpose-of-this-document)
2. [Project Background and Initial Request](#2-project-background-and-initial-request)
3. [User Contributions (Domain Knowledge and Criteria)](#3-user-contributions-domain-knowledge-and-criteria)
4. [AI Agent Contributions](#4-ai-agent-contributions-github-copilot-coding-agent)
5. [Conversation Flow and Collaborative Process](#5-conversation-flow-and-collaborative-process)
6. [Technical Fixes Applied During the Workflow](#6-technical-fixes-applied-during-the-workflow)
7. [Final Outputs and Verification](#7-final-outputs-and-verification)
8. [About the Agent and Model Used](#8-about-the-agent-and-model-used)
9. [How to Run the Code from Scratch (Step-by-Step for Non-Experts)](#9-how-to-run-the-code-from-scratch-step-by-step-for-non-experts)
10. [Data Sources and File Structure](#10-data-sources-and-file-structure)
11. [Glossary](#11-glossary)

---

## 1. Purpose of This Document

This report is intended to make the entire analysis **fully reproducible** by any researcher who has access to:

- ArcGIS Pro with the Spatial Analyst extension, and / or  
- R (≥ 4.1) with the `terra` and `sf` packages

It documents **what was done**, **who contributed what**, **why each decision was made**, and **how to re-run every step from scratch**. A non-expert section (Section 9) provides plain-language instructions for running the workflow without prior GIS programming experience.

---

## 2. Project Background and Initial Request

### 2.1 Research Context

The work belongs to a broader **coral-cover modeling project** for St. Thomas and St. John (USVI), part of the NSF CoPE (Coral Reefs in a Rapidly Changing Environment) program. The goal is to predict benthic coral-cover percentages from environmental predictors, including terrain descriptors derived from high-resolution (2 m) bathymetry data.

Terrain descriptors are used as spatial predictors in a statistical model that relates survey observations of hard-coral cover to measurable environmental variables. The terrain predictors chosen for this project are:

| Predictor | Ecological rationale |
|---|---|
| **Slope (local, 50 m)** | Substrate angle; affects sedimentation and larval settlement |
| **Slope (broad-scale, 240 m)** | Regional exposure to currents and wave energy |
| **Slope-of-slope** | Terrain ruggedness; structural complexity proxy |
| **Aspect** | Directional orientation; related to light and current exposure |
| **Curvature** | Concavity / convexity; controls sediment convergence |

### 2.2 The Initial Question

The user's original request was (paraphrased from the conversation):

> *"I need to build spatial terrain predictors for St. Thomas USVI from a 2 m bathymetry raster. The target coordinate reference system is WGS 84 / UTM Zone 20N (EPSG:32620) and the final output resolution should be 50 m × 50 m. I need slope at multiple scales, slope-of-slope, aspect, and curvature. Can you write a Python script using ArcGIS Pro / arcpy to do this?"*

The user then shared an initial version of the Python script (written using ArcPy) and asked the AI to review, correct, and extend it. Subsequent messages involved iterative debugging as issues were discovered through execution.

---

## 3. User Contributions (Domain Knowledge and Criteria)

The user provided all of the following. Without this domain context, the AI could not have built the correct analysis:

### 3.1 Study Area and Biological Framing
- Study area: **St. Thomas + St. John** (STT/STJ), US Virgin Islands.
- Ecological framing: coral-cover modeling within the NSF CoPE research program.
- List of terrain predictors needed and their ecological justification.

### 3.2 Spatial Specifications
- **Source data:** STTSTJ_2m.tif — a 2 m native-resolution bathymetry raster in the project's Google Drive.
- **Target CRS:** EPSG:32620 (WGS 84 / UTM Zone 20N).
- **Target resolution:** 50 m × 50 m.
- **Focal neighborhood scale for "broad slope":** 240 m radius.
- **Output folder structure** on the shared Google Drive (NSF CoPE internal drive).

### 3.3 Tooling Choice
- The user specified that **ArcGIS Pro + arcpy** should be used (the project team has access to ArcGIS licenses through their institution).
- The user also requested an **R script** for exploratory verification and for the reprojection / resampling step.

### 3.4 Quality-Control Criteria
- All final outputs must be exactly **50.0 m × 50.0 m** cells.
- All final outputs must be in **EPSG:32620** (not EPSG:26920, not geographic CRS).
- Intermediate 2 m files must be kept for QA/archive purposes.
- The file `02_reproject_resample.R` was initially required as a fallback reprojection step; once the Python workflow was corrected to produce true 50 m outputs, the user confirmed that the R script was **no longer needed**.

### 3.5 Execution Feedback (Bug Reports)
The user ran each block interactively in a **Jupyter Notebook** and reported the exact error messages and unexpected output sizes, which allowed the AI to diagnose the root causes of each bug.

---

## 4. AI Agent Contributions (GitHub Copilot Coding Agent)

The AI agent (GitHub Copilot, see Section 8) contributed the following:

### 4.1 Code Architecture
- Designed the block-by-block Jupyter Notebook structure (BLOCK 0 through BLOCK 8) so that each processing step could be run and verified independently without re-running the entire pipeline.
- Created the `WorkflowConfig` class to centralise all configuration parameters (paths, CRS, resolutions, aggregation factors, focal radii) so that changing one value cascades correctly.
- Wrote the `Logger` utility class for timestamped console output.
- Wrote the `_save_raster()` helper to work around silent failures when saving large rasters directly to Google Drive.

### 4.2 Technical Diagnosis and Bug Fixes

Seven critical bugs were identified and corrected (see Section 6 for details). The AI diagnosed:

1. **Wrong aggregation factor** in the original configuration.
2. **Missing CRS reprojection** — terrain derivatives were being computed from the native NAD83 raster instead of the WGS84 target.
3. **Two-step aggregation chain** that left the "50 m" slope raster at 2 m resolution.
4. **Wrong focal radius** — the 240 m neighbourhood was being specified as 5 cells (= 10 m) instead of 120 cells (= 240 m).
5. **Slope-of-slope computed from a mislabelled 2 m raster** instead of the true 50 m raster.
6. **Aspect aggregated with arithmetic mean** — incorrect for a circular variable; replaced with sin/cos circular-mean decomposition.
7. **Curvature not aggregated** — the R reprojection script was being relied upon instead of computing the 50 m output directly in Python.

### 4.3 R Scripts
- Wrote `01_exploratory_spatial_analysis.R`: loads all raster and vector inputs, clips the WCMC coral-reef polygon mask to the study area, reports CRS and resolution of all layers, and counts survey points inside/outside the coral mask.
- Wrote `02_reproject_resample.R`: reprojects and resamples all terrain rasters from EPSG:26920 to EPSG:32620 at 50 m using R's `terra` package, including a circular-mean aspect aggregation. *(This script remains in the repository for reference but is no longer needed now that the Python workflow produces correct 50 m outputs directly.)*

### 4.4 Documentation and Comments
- All code is extensively commented, explaining not just *what* each step does but *why* it is done that way (ecological rationale, statistical correctness, risk of alternative approaches).
- Each BLOCK begins with a docstring that lists the fix applied, the expected outputs, and the estimated processing time.

---

## 5. Conversation Flow and Collaborative Process

### 5.1 How the Interface Works

The work was done using **GitHub Copilot's Coding Agent** (also called the "Task Agent"), which operates directly on a GitHub repository. Here is how the interface works:

1. The user submits a task description (in natural language) via the GitHub Copilot interface.
2. The agent clones the repository into a sandboxed environment, reads all existing files, and plans the work.
3. The agent writes, edits, and commits code changes directly to a feature branch in the repository.
4. The user can see all changes via Pull Requests on GitHub.
5. The agent reports progress at each meaningful step.
6. Conversations are iterative: the user can submit follow-up tasks, error messages, or new requirements, and the agent continues working on the same branch.

### 5.2 Conversation Summary

The following is a high-level summary of the conversation stages:

| Stage | User action | AI action |
|---|---|---|
| **1. Initial request** | Described the research goal, data, and required outputs | Reviewed existing code, identified seven bugs, proposed a corrected architecture |
| **2. First code version** | Ran Block 0 and Block 1; confirmed workspace was set up | Wrote and committed Blocks 0–2 |
| **3. Slope computation** | Ran Block 3; reported unexpected file sizes | Diagnosed wrong aggregation factor; fixed `AGG_FACTOR = 25`, switched from two-step chain to single `Aggregate(factor=25)` |
| **4. 240 m scale slope** | Ran Block 4; reported processing was too fast (seconds instead of minutes) | Diagnosed wrong focal radius (5 cells instead of 120); fixed `FOCAL_RADIUS_240M_at2m = 120` |
| **5. Aspect and curvature** | Ran Blocks 6–7; asked about circular mean | Explained the 0°/360° wrap-around problem for aspect; implemented sin/cos decomposition; added 50 m aggregation for curvature directly in Python |
| **6. Final verification** | Ran Block 8; all five outputs showed ✓ 50.0 m × 50.0 m / EPSG:32620 | Confirmed workflow complete; noted R script is no longer needed |
| **7. Documentation** | Asked for this reproducibility report | Wrote this document |

---

## 6. Technical Fixes Applied During the Workflow

The following table summarises all seven corrections made to the original code. Each fix is labelled in the script header and in the relevant BLOCK docstring.

| Fix | Location | Original (buggy) behaviour | Corrected behaviour |
|---|---|---|---|
| **Fix 1** | `WorkflowConfig` | `AGG_FACTOR` was wrong; two-step factor chain confused cell sizes | Single `AGG_FACTOR = 25` (50 m / 2 m); removed misleading intermediate factors |
| **Fix 2** | BLOCK 1 | Terrain tools ran on native NAD83 bathymetry (EPSG:26920) | Reproject bathymetry to EPSG:32620 first; all derivatives inherit correct CRS |
| **Fix 3** | BLOCK 3 | Two-step chain (Slope → Copy → BlockStatistics) left output at 2 m with wrong filename | Single `Aggregate(factor=25, MEAN)` on slope_2m.tif produces true 50 m cells |
| **Fix 4** | BLOCK 4 | `NbrCircle(5)` on 2 m raster = 10 m neighbourhood | `NbrCircle(120)` on 2 m raster = 240 m neighbourhood |
| **Fix 5** | BLOCK 5 | Slope-of-slope computed from a 2 m raster mislabelled as "50 m" | Uses the true 50 m slope from Fix 3 |
| **Fix 6** | BLOCK 6/6b | Arithmetic mean of aspect degrees — wrong at the 0°/360° wrap | Circular mean: sin/cos decomposition → aggregate → atan2 recovery |
| **Fix 7** | BLOCK 7/7b | Only curvature_2m.tif was produced; aggregation relied on R script | `Aggregate(factor=25, MEAN)` added directly in Python; 50 m output produced here |

---

## 7. Final Outputs and Verification

Block 8 of the Python script performs a final verification and produces the following table (reproduced here from the successful run):

```
File                           CRS / EPSG                          Cell X (m) Cell Y (m)
------------------------------------------------------------------------------------------
  ✓ slope_50m.tif                WGS_1984_UTM_Zone_20N (32620)             50.0       50.0
  ✓ slope_240m_50m.tif           WGS_1984_UTM_Zone_20N (32620)             50.0       50.0
  ✓ slope_of_slope_50m.tif       WGS_1984_UTM_Zone_20N (32620)             50.0       50.0
  ✓ aspect_50m.tif               WGS_1984_UTM_Zone_20N (32620)             50.0       50.0
  ✓ curvature_50m.tif            WGS_1984_UTM_Zone_20N (32620)             50.0       50.0

✓ ALL outputs are in EPSG:32620 at 50m × 50m.
✓ The R script 02_reproject_resample.R is no longer needed.
```

### Intermediate 2 m archives (kept for QA)

| Path (relative to `02_terrain_analysis_outputs/`) | Description |
|---|---|
| `00_bathymetry_source/STTSTJ_2m_native.tif` | Original bathymetry in native CRS |
| `00_bathymetry_source/STTSTJ_2m_EPSG32620.tif` | Bathymetry reprojected to EPSG:32620, 2 m |
| `00_intermediate/slope_2m.tif` | Slope at native 2 m resolution |
| `00_intermediate/slope_2m_240m_focal.tif` | Focal-mean slope at 240 m scale, 2 m pixels |
| `03_aspect/aspect_2m.tif` | Aspect at 2 m resolution (archive) |
| `02_curvature/curvature_2m.tif` | Curvature at 2 m resolution (archive) |

---

## 8. About the Agent and Model Used

| Item | Details |
|---|---|
| **Platform** | GitHub Copilot (github.com) |
| **Agent type** | Copilot Coding Agent / Task Agent |
| **Underlying AI model** | Claude Sonnet (Anthropic) — the version deployed in the GitHub Copilot coding agent at the time of this work (April 2026) |
| **Interface** | Asynchronous task-based chat: the user submits a task description; the agent works autonomously on a branch and reports progress via Pull Requests |
| **Repository integration** | The agent reads, creates, and edits files directly in the GitHub repository; all changes are committed and pushed to a feature branch |
| **Execution environment** | Sandboxed Linux container with internet access; the agent cannot directly run ArcPy (ArcGIS is not available in the sandbox), but it can read, write, and reason about the code |

### Note on AI limitations

The AI agent **did not execute the ArcPy code itself** — it reasoned about correctness from the code, the reported error messages, and the output file sizes that the user provided. All actual raster processing was performed by the user on their local ArcGIS Pro installation. The agent's role was to diagnose bugs, write correct code, and provide explanations.

---

## 9. How to Run the Code from Scratch (Step-by-Step for Non-Experts)

This section explains how to reproduce the full terrain analysis workflow, assuming you have never run a Python or R script before. Follow each numbered step in order.

---

### 9.1 What You Need Before Starting

**Software (all required):**

| Software | Where to get it | Notes |
|---|---|---|
| **ArcGIS Pro** (version 3.x) | Your institution's software portal (requires a license) | Must include the **Spatial Analyst Extension** |
| **Jupyter Notebook** | Installed with ArcGIS Pro or via [conda](https://docs.conda.io/) | Used to run the Python script block by block |
| **R** (version 4.1 or later) | https://cran.r-project.org/ | Free. Only needed for the exploratory R scripts |
| **RStudio** (optional but recommended) | https://posit.co/download/rstudio-desktop/ | Free IDE for R |

**Data (required):**

| File | Description |
|---|---|
| `STTSTJ_2m.tif` | High-resolution (2 m) bathymetry raster for St. Thomas / St. John. Must be placed at the path specified in `WorkflowConfig.BATHYMETRY_SOURCE_TIF` inside the script. |
| WCMC coral-reef shapefile | `WCMC008_CoralReef2021_Py_v4_1.shp` — used only in the R exploratory script |
| Hard coral survey points | `hard_corals_data_filtered.shp` — used only in the R exploratory script |

---

### 9.2 Step 1: Set Up Your Folder Structure

Create the following folders on your computer (or shared drive). The script creates them automatically, but it helps to understand the layout:

```
02_terrain_analysis_outputs/
├── 00_bathymetry_source/    ← input + reprojected bathymetry
├── 00_intermediate/         ← temporary 2m files
├── 01_slope/                ← final slope outputs (50m)
├── 02_curvature/            ← curvature outputs (50m)
└── 03_aspect/               ← aspect outputs (50m)
```

---

### 9.3 Step 2: Open the Python Script in Jupyter Notebook

1. Open **ArcGIS Pro**.
2. In the ribbon at the top, click **Analysis → Python** to open the Python Notebook panel, OR open **Jupyter Notebook** from the Start menu (the ArcGIS Pro conda environment includes Jupyter).
3. Navigate to the folder where you saved `Coral_Modeling_Terrain_Analysis_USVI_202604.py`.
4. Open the file. If it has a `.py` extension, you may need to copy its contents into a new Jupyter Notebook (`.ipynb`) manually — create one cell per BLOCK section.

> **Tip:** If you are unsure how to open Jupyter Notebook inside ArcGIS Pro, search for "ArcGIS Pro Python Notebook" in your ArcGIS Pro help documentation.

---

### 9.4 Step 3: Edit the File Paths

Before running any block, open the script and **update the two path variables** to match your system:

```python
# In WorkflowConfig class (top of the script):
PROJECT_ROOT = r"G:\Shared drives\NSF CoPE internal\..."   # ← change this to YOUR path
BATHYMETRY_SOURCE_TIF = r"G:\...\STTSTJ_2m.tif"           # ← change this to YOUR path
```

Use the **full path** to the folder and file on your computer. On Windows, paths use backslashes (`\`). The leading `r` before the quote means "raw string" — keep it there so that backslashes are interpreted correctly.

---

### 9.5 Step 4: Run the Blocks in Order

Run each BLOCK cell one at a time, in order. **Do not skip blocks** on the first run. Wait for each block to finish before running the next one.

| BLOCK | What it does | Expected time |
|---|---|---|
| **BLOCK 0** | Loads configuration; checks that ArcGIS extensions are available | < 5 seconds |
| **BLOCK 1** | Creates output folders; reprojects bathymetry to EPSG:32620 | 1–2 minutes |
| **BLOCK 2** | Validates that the bathymetry file exists and is readable | < 5 seconds |
| **Verification 2A** | Prints statistics about the bathymetry (depth range, CRS, cell size) | < 5 seconds |
| **BLOCK 3** | Computes slope at 2 m, then aggregates to true 50 m | 1–2 minutes |
| **BLOCK 4** | Computes slope at 240 m scale (focal mean) — this is slow | **5–8 minutes** |
| **BLOCK 5** | Computes slope-of-slope from the true 50 m slope | 30–40 seconds |
| **BLOCK 6 + 6b** | Computes aspect at 2 m and aggregates to 50 m using circular mean | 2–3 minutes |
| **BLOCK 7 + 7b** | Computes curvature at 2 m and aggregates to 50 m | 1–2 minutes |
| **BLOCK 8** | Verifies that all five final outputs are correct (50 m, EPSG:32620) | < 10 seconds |

Each block prints messages beginning with `✓` (success) or `⚠` (warning) to the console. **Do not proceed to the next block if you see a `✗` (error) message.**

---

### 9.6 Step 5: Verify the Final Outputs

After BLOCK 8 runs successfully, you should see exactly this output:

```
  ✓ slope_50m.tif                WGS_1984_UTM_Zone_20N (32620)   50.0   50.0
  ✓ slope_240m_50m.tif           WGS_1984_UTM_Zone_20N (32620)   50.0   50.0
  ✓ slope_of_slope_50m.tif       WGS_1984_UTM_Zone_20N (32620)   50.0   50.0
  ✓ aspect_50m.tif               WGS_1984_UTM_Zone_20N (32620)   50.0   50.0
  ✓ curvature_50m.tif            WGS_1984_UTM_Zone_20N (32620)   50.0   50.0
```

All five files are now ready to use as spatial predictors in your coral-cover model.

---

### 9.7 Step 6 (Optional): Run the R Exploratory Script

If you want to run the exploratory spatial analysis in R:

1. Open **RStudio**.
2. Open `01_exploratory_spatial_analysis.R`.
3. Edit the `root` variable at the top to match your folder path.
4. Install required packages by running this line in the RStudio console:
   ```r
   install.packages(c("terra", "sf", "dplyr"))
   ```
5. Run the entire script (Ctrl+Shift+Enter in RStudio, or click **Run** → **Run All**).

The script will:
- Clip the global WCMC coral-reef polygon to the study area.
- Report the CRS and resolution of all raster and vector files.
- Count how many hard-coral survey points fall inside the coral-reef mask.

---

### 9.8 Troubleshooting Common Errors

| Error message | Likely cause | Solution |
|---|---|---|
| `FileNotFoundError: STTSTJ_2m.tif not found` | The path in `WorkflowConfig` is wrong | Double-check the path and make sure the file exists |
| `RuntimeError: Spatial Analyst Extension unavailable` | Extension is not licensed or checked out | Open ArcGIS Pro → Project → Licensing → Check out Spatial Analyst |
| `PermissionError` writing to Google Drive | Google Drive sync conflict | Save to a local folder first, then copy to Drive |
| BLOCK 4 finishes in < 10 seconds | The focal radius is wrong (should be 120 cells) | Make sure `FOCAL_RADIUS_240M_at2m = 120` in `WorkflowConfig` |
| Output cell size is 2.0 m instead of 50.0 m | Wrong `AGG_FACTOR` | Make sure `AGG_FACTOR = 25` in `WorkflowConfig` |

---

## 10. Data Sources and File Structure

### 10.1 Primary Input Data

| Dataset | Source | Native CRS | Resolution |
|---|---|---|---|
| STTSTJ_2m.tif — bathymetry | NOAA / NSF CoPE project | NAD83 / UTM Zone 20N (EPSG:26920) or WGS84/UTM 20N | 2 m |
| WCMC008_CoralReef2021_Py_v4_1.shp | UNEP-WCMC 2021 coral atlas | Geographic WGS84 | vector |
| hard_corals_data_filtered.shp | NSF CoPE field surveys, St. Thomas | Varies (reprojected in script) | vector |

### 10.2 Output File Naming Convention

| Pattern | Meaning |
|---|---|
| `*_2m.tif` | Computed at native 2 m resolution (intermediate or archive) |
| `*_50m.tif` | Final product; every cell is a true 50 m × 50 m pixel in EPSG:32620 |
| `*_240m_50m.tif` | Slope smoothed over a 240 m neighbourhood, aggregated to 50 m cells |

---

## 11. Glossary

| Term | Definition |
|---|---|
| **EPSG code** | A numerical identifier for a coordinate reference system. EPSG:32620 = WGS 84 / UTM Zone 20N (used for most of the Caribbean) |
| **CRS (Coordinate Reference System)** | The mathematical framework that defines how map coordinates relate to locations on the Earth |
| **UTM** | Universal Transverse Mercator — a projected CRS that uses metres as units, making distance calculations straightforward |
| **Resolution / cell size** | The side length of each pixel in a raster. 50 m means each pixel covers a 50 m × 50 m square on the ground |
| **Aggregation** | Combining multiple small pixels into one larger pixel (e.g., 25×25 = 625 pixels at 2 m merged into one pixel at 50 m) |
| **Focal statistics** | A moving-window operation where each pixel is replaced by a summary statistic (e.g., mean) computed over a neighbourhood of cells |
| **Aspect** | The compass direction (0°–360°) that a slope face points toward |
| **Curvature** | The rate of change of slope. Positive = convex (ridgeline); Negative = concave (valley) |
| **Circular mean** | The correct way to average angular measurements (like aspect) that wrap around at 0°/360° |
| **arcpy** | The Python package bundled with ArcGIS Pro that provides access to all ArcGIS geoprocessing tools |
| **terra** | An R package for raster and vector spatial analysis |
| **Snap raster** | A reference raster that forces all outputs to align to the same grid origin and cell boundaries |

---

*This document was generated by the GitHub Copilot Coding Agent (Claude Sonnet) in collaboration with J. Olaya, April 2026.*  
*It should be updated if the code or workflow changes significantly.*
