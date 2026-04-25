# Coral Cover Restoration and Bleaching Scenario Modeling
# Spatial Modeling Workflow – Belize Coral Reef System
###############################################################

# This script models coral reef condition under sequential restoration and bleaching scenarios using raster data

# -------------------------------------------------------------
# Load required libraries
# -------------------------------------------------------------
library(terra)       # modern spatial raster package
library(tidyverse)   # for data manipulation + plotting
# -------------------------------------------------------------

# -------------------------------------------------------------
# Setup Directories
# -------------------------------------------------------------
# Set your working directory to the project root or use relative paths
# Suggested structure: /your_project_folder/data/ and /your_project_folder/outputs/
input_dir  <- "path/to/your/input_data"
output_dir <- "path/to/your/outputs"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------------------------------------------
# Load the Spatial Data
# -------------------------------------------------------------
# Ensure the rasters and shapefiles are placed in your data directory
baseline      <- rast(file.path(input_dir, "cor_baseline.tif"))
depth         <- rast(file.path(input_dir, "depth.tif"))
fishing_zones <- vect(file.path(input_dir, "ma2016.shp"))
fishing_zones_utm <- project(fishing_zones, crs(baseline))

# Inspect metadata
baseline
depth
fishing_zones

# =============================================================
# Align depth raster with coral raster and build shallow mask
# =============================================================
#  Crop depth raster to baseline extent
depth_crop <- crop(depth, baseline)

#  Resample depth to match baseline grid
depth_res <- resample(depth_crop, baseline)

writeRaster(depth_res, file.path(out_dir, "depth_resampled_baseline.tif"), overwrite=TRUE)

# =============================================================
# Create Depth Masks
# logical rasters are later used to apply bleaching intensity
# =============================================================
# shallow reef 0–5 m (>= -5 m)
mask_0_5m  <- depth_res >= -5

# mid-depth reef 5–10 m (>= -10 m AND < -5 m)
mask_5_10m <- depth_res >= -10 & depth_res < -5

# deep reef >10 m (< -10 m); bleaching attenuated ~12% loss (Hughes et al. 2018, doi:10.3354/meps12732)
mask_10_27m <- depth_res < -10

# sanity: resample to coral grid
mask_0_5m    <- resample(mask_0_5m,   baseline, method="near")
mask_5_10m   <- resample(mask_5_10m,  baseline, method="near")
mask_10_27m  <- resample(mask_10_27m, baseline, method="near")

#---- define fishing zone names
zone_field <- "Name"

# -------------------------------------------------------------
# Function definition
#--------------------------------------------------------------

# -----------
# Function to classify degraded (<15%) and healthy (≥15%) coral
classify_coral <- function(r, threshold = 15) {
  classify(
    r,
    rcl = matrix(c(
      -Inf, threshold, 1,   # 1 = degraded
      threshold, Inf, 2     # 2 = healthy
    ), ncol = 3, byrow = TRUE)
  )
}

#----------
# Function to Compute area and % for a categorical raster (1=degraded, 2=healthy)
# Pixel size = 10 m × 10 m = 100 m² = 0.0001 km²
area_by_class <- function(r_class) {
  
  pixel_area_m2 <- prod(res(r_class))   # 10 × 10 = 100 m²
  tab <- freq(r_class)                  # terra version, no useNA argument
  
  total_pixels <- sum(tab$count)
  
  tibble(
    class = c("Degraded (<15%)", "Healthy (≥15%)"),
    pixels = tab$count,
    percent = (tab$count / total_pixels) * 100,
    area_m2 = tab$count * pixel_area_m2,
    area_km2 = (tab$count * pixel_area_m2) / 1e6
  )
}

#----------
# Apply restoration: increase coral cover by +10% in degraded pixels
# Relative increase in coral cover ONLY for degraded pixels
# factor = 0.10 means +10% relative; factor = 0.25 means +25% relative
increase_coral_relative <- function(r, factor, degraded_mask) {
  
  # Use terra::ifel for pixel-wise conditional:
  # if degraded (mask == 1) → r * (1 + factor)
  # else                  → r
  r_new <- ifel(degraded_mask == 1, r * (1 + factor), r)
  
  # Cap values at 100%
  r_new[r_new > 100] <- 100
  
  return(r_new)
}

#----------
# ============================================================
# NEW bleaching function:
#   0–5 m   : 50% loss
#   5–10 m  : 25% loss
#   >10 m   : 12% loss  (Hughes et al. 2018, doi:10.3354/meps12732)
# ============================================================

apply_bleaching_depth <- function(r, mask_0_5m, mask_5_10m, mask_10_27m){
  
  # Start copy
  r_new <- r
  
  # 50% loss for 0–5 m
  r_new[mask_0_5m] <- r[mask_0_5m] * (1 - 0.50)
  
  # 25% loss for 5–10 m
  r_new[mask_5_10m] <- r[mask_5_10m] * (1 - 0.25)
  
  # 12% loss for >10 m (curvilinear attenuation, deep mesophotic)
  r_new[mask_10_27m] <- r[mask_10_27m] * (1 - 0.12)
  
  # No negatives
  r_new[r_new < 0] <- 0
  
  return(r_new)
}


#-----------------------------------------------------------
# baseline description
#-----------------------------------------------------------
# Summarize coral condition by fishing zone
# For each polygon in fishing_zones, compute area degraded/healthy within that zone.
summarize_coral_by_zones <- function(r_class, zones, zone_field, scenario_name) {
  
  pixel_area_m2 <- prod(res(r_class))
  
  # Loop over each fishing zone polygon
  zone_results <- lapply(1:nrow(zones), function(i) {
    
    this_zone <- zones[i, ]
    zone_id   <- as.character(this_zone[[zone_field]])
    
    # Mask raster to the current fishing zone (cells outside become NA)
    r_zone <- mask(r_class, this_zone)
    
    tab <- freq(r_zone)
    
    # If no coral pixels in this zone, skip (returns NULL)
    if (is.null(tab)) return(NULL)
    
    total_pixels <- sum(tab$count)
    
    tibble(
      scenario     = scenario_name,
      zone         = zone_id,
      class        = c("Degraded (<15%)", "Healthy (≥15%)"),
      pixels       = tab$count,
      percent_zone = (tab$count / total_pixels) * 100,  # % within this zone
      area_m2      = tab$count * pixel_area_m2,
      area_km2     = (tab$count * pixel_area_m2) / 1e6
    )
  })
  
  bind_rows(zone_results)
}

# =============================================================
# Baseline description (national level)
# =============================================================

baseline_class <- classify_coral(baseline, threshold = 15)

writeRaster(baseline_class, file.path(out_dir, "baseline_class.tif"), overwrite=TRUE)

# Baseline map
plot(
  baseline_class,
  main = "Baseline Coral Condition (1 = Degraded, 2 = Healthy)",
  col = c("red", "turquoise3"),
  axes = TRUE
)

# National baseline area
baseline_area <- area_by_class(baseline_class)
baseline_area

# Baseline barplot (national)
ggplot(baseline_area,
       aes(x = class, y = area_km2, fill = class)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Degraded (<15%)" = "red",
    "Healthy (≥15%)"  = "turquoise3"
  )) +
  labs(title = "Area of Degraded vs Healthy Coral (Baseline)",
       x = "Coral Condition",
       y = "Area (km²)",
       fill = "Coral Condition") +
  theme_minimal()


# =============================================================
# Bleaching Only Scenario (no restoration)
# Applies bleaching event directly on the baseline raster
# (no Phase 1 or Phase 2 restoration)
# =============================================================

bleaching_only <- apply_bleaching_depth(
  r = baseline,
  mask_0_5m   = mask_0_5m,
  mask_5_10m  = mask_5_10m,
  mask_10_27m = mask_10_27m
)

writeRaster(bleaching_only, file.path(out_dir, "bleaching_only_coral_cover.tif"), overwrite=TRUE)

# Reclassify Bleaching Only
bleaching_only_class <- classify_coral(bleaching_only, threshold = 15)

writeRaster(bleaching_only_class, file.path(out_dir, "bleaching_only_class.tif"), overwrite=TRUE)

plot(
  bleaching_only_class,
  col = c("red", "turquoise3"),
  main = "Bleaching Only Coral Condition (1 = Degraded, 2 = Healthy)",
  axes = TRUE
)

# National area for Bleaching Only
bleaching_only_area <- area_by_class(bleaching_only_class)
bleaching_only_area

# Plot Bleaching Only areas
ggplot(bleaching_only_area,
       aes(x = class, y = area_km2, fill = class)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Degraded (<15%)" = "red",
    "Healthy (\u226515%)"  = "turquoise3"
  )) +
  labs(title = "Area of Degraded vs Healthy Coral (Bleaching Only)",
       x = "Coral Condition",
       y = "Area (km\u00b2)",
       fill = "Coral Condition") +
  theme_minimal()

# Loss map: Baseline healthy pixels that become degraded after bleaching only
bleaching_only_loss <- baseline_class == 2 & bleaching_only_class == 1

writeRaster(bleaching_only_loss, file.path(out_dir, "bleaching_only_loss.tif"), overwrite=TRUE)

plot(
  bleaching_only_loss,
  col = c("lightgrey", "red"),
  main = "Loss of Healthy Coral \u2014 Bleaching Only Scenario",
  axes = TRUE
)

# =====================================================================
# Phase 1 — Short-term restoration (+10% coral cover on degraded corals)
# =====================================================================

# Mask: degraded in baseline = 1; healthy = 2
baseline_degraded_mask <- baseline_class == 1

# Apply Phase 1 restoration
phase1 <- increase_coral_relative(
  r = baseline,
  factor = 0.10,   # +10% relative increase
  degraded_mask = baseline_degraded_mask
)

writeRaster(phase1, file.path(out_dir,"phase1_coral_cover.tif"), overwrite=TRUE)


# Reclassify Phase 1
phase1_class <- classify_coral(phase1, threshold = 15)

writeRaster(phase1_class, file.path(out_dir,"phase1_class.tif"), overwrite=TRUE)

plot(
  phase1_class,
  col = c("red", "turquoise3"),
  main = "Phase 1 Coral Condition (1 = Degraded, 2 = Healthy)",
  axes = TRUE
)

# National area after Phase 1
phase1_area <- area_by_class(phase1_class)
phase1_area

# Plot Phase 1 areas
ggplot(phase1_area,
       aes(x = class, y = area_km2, fill = class)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Degraded (<15%)" = "red",
    "Healthy (≥15%)"  = "turquoise3"
  )) +
  labs(title = "Area of Degraded vs Healthy Coral (Phase 1)",
       x = "Coral Condition",
       y = "Area (km²)",
       fill = "Coral Condition") +
  theme_minimal()

# Map of improvements: baseline degraded → Phase 1 healthy
phase1_improved <- baseline_class == 1 & phase1_class == 2

writeRaster(phase1_improved, file.path(out_dir, "phase1_improved.tif"), overwrite=TRUE)

plot(
  phase1_improved,
  col = c("lightgrey", "turquoise3"),
  main = "Pixels Improved After Phase 1 (Degraded → Healthy)",
  axes = TRUE
)

# =============================================================
# First bleaching event (on Phase 1)
# =============================================================
phase1_bleach <- apply_bleaching_depth(
  r = phase1,
  mask_0_5m   = mask_0_5m,
  mask_5_10m  = mask_5_10m,
  mask_10_27m = mask_10_27m
)

writeRaster(phase1_bleach, file.path(out_dir, "phase1_bleach_coral_cover.tif"), overwrite=TRUE)

# Reclassify Phase 1 after bleaching
phase1_bleach_class <- classify_coral(phase1_bleach, 15)

writeRaster(phase1_bleach_class, file.path(out_dir, "phase1_bleach_class.tif"), overwrite=TRUE)

plot(
  phase1_bleach_class,
  col = c("red", "turquoise3"),
  main = "Phase 1 After Bleaching (1 = Degraded, 2 = Healthy)",
  axes = TRUE
)

# National area after Phase 1 bleaching
phase1_bleach_area  <- area_by_class(phase1_bleach_class)
phase1_bleach_area

# Loss map: Healthy → Degraded (bleaching loss)
phase1_loss <- phase1_class == 2 & phase1_bleach_class == 1

writeRaster(phase1_loss, file.path(out_dir, "phase1_loss.tif"), overwrite=TRUE)

plot(
  phase1_loss,
  col = c("lightgrey", "red"),
  main = "Loss of Healthy Coral After Bleaching (Phase 1)",
  axes = TRUE
)

# =============================================================
# Phase 2 — Medium-term restoration (+25% relative on Phase1-bleached)
# Restoration is applied ONLY where baseline was degraded.
# =============================================================

# first without climate change (increased per phase1)
phase2 <- increase_coral_relative(
  r = phase1_bleach,           # starting from Phase 1 AFTER bleaching
  factor = 0.25,               # +25% relative increase
  degraded_mask = baseline_degraded_mask
)

phase2[phase2 > 100] <- 100    # safety cap at 100%

writeRaster(phase2, file.path(out_dir, "phase2_coral_cover.tif"), overwrite=TRUE)

# Reclassify Phase 2
phase2_class <- classify_coral(phase2, threshold = 15)

writeRaster(phase2_class, file.path(out_dir, "phase2_class.tif"), overwrite=TRUE)

plot(
  phase2_class,
  col = c("red", "turquoise3"),
  main = "Phase 2 Coral condition (1 = Degraded, 2 = Healthy)",
  axes = TRUE
)

# Area
phase2_area <- area_by_class(phase2_class)
phase2_area

# Improvements map
phase2_improved <- baseline_class == 1 & phase2_class == 2

writeRaster(phase2_improved, file.path(out_dir, "phase2_improved.tif"), overwrite=TRUE)

plot(
  phase2_improved,
  col = c("lightgrey", "turquoise3"),
  main = "Pixels Improved After Phase 2 (Degraded → Healthy)",
  axes = TRUE
)

### area estimation
library(terra)

phase2_improved
plot(phase2_improved, main = "Phase 2 Improved Pixels")

zones_all <- bind_rows(
  baseline_zones,
  phase1_zones,
  phase1_bleach_zones,
  phase2_zones,
  phase2_bleach_zones
)

summarize_phase2_gains_by_zone <- function(change_raster, zones, zone_field) {
  
  pixel_area_m2 <- prod(res(change_raster))
  
  results <- lapply(1:nrow(zones), function(i) {
    
    this_zone <- zones[i, ]
    zone_id   <- as.character(this_zone[[zone_field]])
    
    # Mask change raster to zone
    change_zone <- mask(change_raster, this_zone)
    
    # Count improved pixels (TRUE values)
    improved_pixels <- global(change_zone, "sum", na.rm = TRUE)[1,1]
    
    if (is.na(improved_pixels)) improved_pixels <- 0
    
    tibble(
      zone = zone_id,
      improved_pixels = improved_pixels,
      improved_area_m2 = improved_pixels * pixel_area_m2,
      improved_area_km2 = (improved_pixels * pixel_area_m2) / 1e6
    )
  })
  
  bind_rows(results)
}

phase2_gains_by_zone <- summarize_phase2_gains_by_zone(
  change_raster = phase2_improved,
  zones = fishing_zones_utm,
  zone_field = zone_field
)

phase2_gains_by_zone

# Percentage of Zone Improved
# Get baseline degraded pixels per zone
baseline_degraded_by_zone <- zones_all |>
  filter(scenario == "Baseline", class == "Degraded") |>
  select(zone, baseline_degraded_pixels = pixels)

# Join
phase2_gains_by_zone <- phase2_gains_by_zone |>
  left_join(baseline_degraded_by_zone, by = "zone") |>
  mutate(
    percent_recovered_of_degraded =
      (improved_pixels / baseline_degraded_pixels) * 100
  )

phase2_gains_by_zone

#####
library(ggplot2)

ggplot(phase2_gains_by_zone,
       aes(x = reorder(zone, improved_area_km2),
           y = improved_area_km2)) +
  geom_col(fill = "turquoise3", color = "black", width = 0.7) +
  coord_flip() +
  labs(
    title = "Recovered Coral Area by Fishing Zone",
    subtitle = "Pixels Degraded in Baseline and Healthy After Phase 2",
    x = "Fishing Zone",
    y = "Recovered Area (km²)"
  ) +
  theme_minimal(base_size = 14)

# =============================================================
# Bleaching Only losses by zone
# Pixels that were Healthy in Baseline and became Degraded
# after the Bleaching Only scenario
# =============================================================

bleaching_only_loss_raster <- rast(
  file.path(out_dir, "spatial_files/bleaching_only_loss.tif")
)

summarize_losses_by_zone <- function(change_raster, zones, zone_field) {

  pixel_area_m2 <- prod(res(change_raster))

  results <- lapply(1:nrow(zones), function(i) {

    this_zone <- zones[i, ]
    zone_id   <- as.character(this_zone[[zone_field]])

    change_zone   <- mask(change_raster, this_zone)
    lost_pixels   <- global(change_zone, "sum", na.rm = TRUE)[1, 1]

    if (is.na(lost_pixels)) lost_pixels <- 0

    tibble(
      zone            = zone_id,
      lost_pixels     = lost_pixels,
      lost_area_m2    = lost_pixels * pixel_area_m2,
      lost_area_km2   = (lost_pixels * pixel_area_m2) / 1e6
    )
  })

  bind_rows(results)
}

bleaching_only_losses_by_zone <- summarize_losses_by_zone(
  change_raster = bleaching_only_loss_raster,
  zones         = fishing_zones_utm,
  zone_field    = zone_field
)

bleaching_only_losses_by_zone

# =============================================================
# Diverging bar plot: Bleaching losses (left, red) vs
# Phase 2 restoration gains (right, blue) by fishing zone
# =============================================================

diverging_data <- bind_rows(
  phase2_gains_by_zone |>
    select(zone, area_km2 = improved_area_km2) |>
    mutate(direction = "Restored (Phase 2)",
           bar_value  = area_km2),

  bleaching_only_losses_by_zone |>
    select(zone, area_km2 = lost_area_km2) |>
    mutate(direction = "Lost (Bleaching Only)",
           bar_value  = -area_km2)
)

# Order zones by net change (gains – losses)
zone_order <- diverging_data |>
  group_by(zone) |>
  summarise(net = sum(bar_value), .groups = "drop") |>
  arrange(net) |>
  pull(zone)

diverging_data <- diverging_data |>
  mutate(zone = factor(zone, levels = zone_order))

ggplot(diverging_data,
       aes(x = zone, y = bar_value, fill = direction)) +
  geom_col(color = "black", width = 0.7) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "Lost (Bleaching Only)"  = "#D62728",
      "Restored (Phase 2)"     = "turquoise3"
    )
  ) +
  coord_flip() +
  labs(
    title    = "Coral Area Change by Fishing Zone",
    subtitle = "Red: lost under Bleaching Only · Blue: recovered under Phase 2 restoration",
    x        = "Fishing Zone",
    y        = "Area (km²)  ←  Lost  |  Recovered  →",
    fill     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

# =============================================================
# Second bleaching event (on Phase 2)
# =============================================================

# Bleaching on Phase 2
phase2_bleach <- apply_bleaching_depth(
  r = phase2,
  mask_0_5m   = mask_0_5m,
  mask_5_10m  = mask_5_10m,
  mask_10_27m = mask_10_27m
)

writeRaster(phase2_bleach, file.path(out_dir, "phase2_bleach_coral_cover.tif"), overwrite=TRUE)

# classify after bleaching event
phase2_bleach_class <- classify_coral(phase2_bleach, 15)

writeRaster(phase2_bleach_class, file.path(out_dir, "phase2_bleach_class.tif"), overwrite=TRUE)

plot(
  phase2_bleach_class,
  col = c("red", "turquoise3"),
  main = "Phase 2 after bleaching event (1 = Degraded, 2 = Healthy)",
  axes = TRUE
)

phase2_bleach_area  <- area_by_class(phase2_bleach_class)
phase2_bleach_area

# Loss map: Healthy → Degraded (Phase 2 bleaching)
phase2_loss <- phase2_class == 2 & phase2_bleach_class == 1

writeRaster(phase2_loss, file.path(out_dir, "phase2_loss.tif"), overwrite=TRUE)

plot(
  phase2_loss,
  col = c("lightgrey", "red"),
  main = "Loss of Healthy Coral After Bleaching (Phase 2)",
  axes = TRUE
)

# =============================================================
# National stacked barplots across scenarios (already done)
# =============================================================

scenario_area_long <- bind_rows(
  baseline_area       |> mutate(scenario = "Baseline"),
  bleaching_only_area |> mutate(scenario = "Bleaching Only"),
  phase1_area         |> mutate(scenario = "Phase 1"),
  phase1_bleach_area  |> mutate(scenario = "Bleaching"),
  phase2_area         |> mutate(scenario = "Phase 2")
) |>
  select(scenario, class, area_km2, percent) |>
  mutate(
    class = recode(class,
                   "Degraded (<15%)" = "Degraded",
                   "Healthy (≥15%)"  = "Healthy"),
    scenario = factor(
      scenario,
      levels = c("Baseline", "Bleaching Only", "Phase 1", "Bleaching", "Phase 2")
    )
  )

library(ragg)
library(scales)

# Build two-panel data for national stacked barplot
scenario_area_2panel <- bind_rows(
  scenario_area_long %>%
    filter(scenario %in% c("Baseline", "Bleaching Only")) %>%
    mutate(scenario_group = "Bleaching Only"),
  scenario_area_long %>%
    filter(scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(scenario_group = "With Restoration")
) %>%
  mutate(scenario_group = factor(scenario_group, levels = c("Bleaching Only", "With Restoration")))

p_2panel <- ggplot(scenario_area_2panel, aes(x = scenario, y = percent, fill = class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Degraded" = "red", "Healthy" = "turquoise3")) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  labs(
    x = "Scenario",
    y = "Percent of Total Coral Reef Area",
    fill = "Coral Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    # --- Legend arriba para optimizar espacio ---
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),

    # reduce un poco el espacio extra arriba
    plot.margin = margin(t = 6, r = 8, b = 6, l = 8),
    strip.text = element_text(face = "bold", size = 13)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p_2panel

ggsave(
  filename = file.path(out_dir, "coral_cover_scenarios_stacked_barplot.png"),
  plot = p_2panel,
  device = ragg::agg_png,
  width = 14, height = 7, units = "in",
  dpi = 300,
  bg = "white"
)

## Two-panel national line plot: Bleaching Only vs With Restoration
library(dplyr)
library(ggplot2)
library(scales)

# 1) Prepare two-panel healthy area data (national level)
healthy_area_2panel <- bind_rows(
  scenario_area_long %>%
    filter(class == "Healthy", scenario %in% c("Baseline", "Bleaching Only")) %>%
    mutate(scenario_group = "Bleaching Only"),
  scenario_area_long %>%
    filter(class == "Healthy", scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(scenario_group = "With Restoration")
) %>%
  mutate(
    scenario_group = factor(scenario_group, levels = c("Bleaching Only", "With Restoration")),
    scenario = factor(scenario, levels = c("Baseline", "Bleaching Only", "Phase 1", "Bleaching", "Phase 2"))
  ) %>%
  arrange(scenario_group, scenario)

# 2) Plot
p_line_2panel <- ggplot(healthy_area_2panel,
                         aes(x = scenario, y = area_km2, group = scenario_group)) +
  geom_line(linewidth = 1.1, color = "#1B9E77") +
  geom_point(size = 3.2, color = "#1B9E77") +
  geom_text(
    aes(label = number(area_km2, accuracy = 1)),
    vjust = -0.8,
    size = 3.6,
    color = "#2B2B2B"
  ) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    title = "Total area of Healthy coral reef (coral cover >15%)",
    subtitle = "Total area (km²)",
    x = NULL,
    y = expression(paste("Healthy coral area (", km^2, ")")),
    caption = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12.5, color = "grey30"),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(color = "grey40", size = 10),
    strip.text = element_text(face = "bold", size = 13)
  )

p_line_2panel

out_dir_fig <- "set/you/output/file"

library(svglite)
ggsave(
   filename = file.path(out_dir_fig, "Fig_HealthyCoralArea_byScenario_Belize.svg"),
   plot = p_line_2panel,
   device = svglite::svglite,
   width = 12, height = 5.2, units = "in"
 )

# =============================================================
# Summarize coral condition by fishing zones (corrected CRS)
# =============================================================

baseline_zones <- summarize_coral_by_zones(
  baseline_class,
  fishing_zones_utm,   # << Correct CRS
  zone_field,
  "Baseline"
)

phase1_zones <- summarize_coral_by_zones(
  phase1_class,
  fishing_zones_utm,   # << Correct CRS
  zone_field,
  "Phase 1"
)

phase1_bleach_zones <- summarize_coral_by_zones(
  phase1_bleach_class,
  fishing_zones_utm,   # << Correct CRS
  zone_field,
  "Bleaching"
)

phase2_zones <- summarize_coral_by_zones(
  phase2_class,
  fishing_zones_utm,   # << Correct CRS
  zone_field,
  "Phase 2"
)

phase2_bleach_zones <- summarize_coral_by_zones(
  phase2_bleach_class,
  fishing_zones_utm,   # << Correct CRS
  zone_field,
  "Phase 2 Bleach"
)

bleaching_only_zones <- summarize_coral_by_zones(
  bleaching_only_class,
  fishing_zones_utm,   # << Correct CRS
  zone_field,
  "Bleaching Only"
)


# Build the combined table
zones_all <- bind_rows(
  baseline_zones,
  bleaching_only_zones,
  phase1_zones,
  phase1_bleach_zones,
  phase2_zones,
  phase2_bleach_zones
) |>
  mutate(
    class = recode(class,
                   "Degraded (<15%)" = "Degraded",
                   "Healthy (≥15%)"  = "Healthy"),
    scenario = recode(as.character(scenario), "Phase 1 Bleach" = "Bleaching"),
    scenario = factor(
      scenario,
      levels = c("Baseline", "Bleaching Only", "Phase 1", "Bleaching",
                 "Phase 2", "Phase 2 Bleach")
    )
  )

zones_all
write.csv(zones_all, file = "your/folder/zones_all_restoration_outputs.csv")

#-------------------- Plots
# Two-panel line plot: Healthy coral (%) across scenarios (by zone)
# Panel 1: Bleaching Only scenario  |  Panel 2: With Restoration pathway
healthy_trend <- zones_all |>
  filter(class == "Healthy", scenario != "Phase 2 Bleach") |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("Baseline", "Bleaching Only", "Phase 1", "Bleaching", "Phase 2")
    )
  )

healthy_trend_2panel <- bind_rows(
  healthy_trend |>
    filter(scenario %in% c("Baseline", "Bleaching Only")) |>
    mutate(scenario_group = "Bleaching Only"),
  healthy_trend |>
    filter(scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) |>
    mutate(scenario_group = "With Restoration")
) |>
  mutate(scenario_group = factor(scenario_group, levels = c("Bleaching Only", "With Restoration")))

p_zone_pct_2panel <- ggplot(healthy_trend_2panel,
       aes(x = scenario, y = percent_zone, group = zone, color = zone)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Trajectories of Healthy Coral (%) by Fishing Zone",
    subtitle = "Effects of restoration and bleaching over scenarios",
    x = NULL,
    y = "Healthy Coral (%)",
    color = "Fishing Zone"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13)
  )

p_zone_pct_2panel

# Healthy coral area by Fishing Zone across scenarios (excluding "Phase 2 Bleach")
zones_healthy <- zones_all %>%
  filter(class == "Healthy", scenario != "Phase 2 Bleach") %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("Baseline", "Bleaching Only", "Phase 1", "Bleaching", "Phase 2")
    ),
    zone = factor(zone)  # keeps a stable legend order
  ) %>%
  arrange(zone, scenario)

### estimate national area and inlcuded in the plot
library(dplyr)
library(ggplot2)
library(scales)
library(forcats) # Required for handling factor levels

### Build two-panel data for zone plots with temporal x-axis labels
# Bleaching Only panel: 4 time points where
#   "1 Year"   = Baseline values (no change in this phase)
#   "2 Years"  = Bleaching Only values (bleaching impact)
#   "~5 Years" = Bleaching Only values (assumed no further change)
# With Restoration panel: same values, renamed to temporal labels
zones_healthy_2panel <- bind_rows(

  # ---- Bleaching Only panel (4 synthetic time points) ----
  zones_healthy %>%
    filter(scenario == "Baseline") %>%
    mutate(scenario = "Baseline", scenario_group = "Bleaching Only"),
  zones_healthy %>%
    filter(scenario == "Baseline") %>%
    mutate(scenario = "1 Year", scenario_group = "Bleaching Only"),
  zones_healthy %>%
    filter(scenario == "Bleaching Only") %>%
    mutate(scenario = "2 Years", scenario_group = "Bleaching Only"),
  zones_healthy %>%
    filter(scenario == "Bleaching Only") %>%
    mutate(scenario = "~5 Years", scenario_group = "Bleaching Only"),

  # ---- With Restoration panel (renamed to temporal labels) ----
  zones_healthy %>%
    filter(scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(
      scenario = recode(as.character(scenario),
                        "Baseline"  = "Baseline",
                        "Phase 1"   = "1 Year",
                        "Bleaching" = "2 Years",
                        "Phase 2"   = "~5 Years"),
      scenario_group = "With Restoration"
    )

) %>%
  mutate(
    scenario       = factor(scenario, levels = c("Baseline", "1 Year", "2 Years", "~5 Years")),
    scenario_group = factor(scenario_group, levels = c("Bleaching Only", "With Restoration"))
  )

### National area per scenario and group (for dashed overlay line)
national_summary_2panel <- zones_healthy_2panel %>%
  group_by(scenario_group, scenario) %>%
  summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
  mutate(zone = "National")

### Zone area line plot (two panels)
p_zone_area_2panel <- ggplot(zones_healthy_2panel,
                              aes(x = scenario, y = area_km2, group = zone, color = zone)) +
  # Regional layers (Individual Fishing Zones)
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +

  # National scale layer (Dashed black line) — linetype mapped for legend entry
  # 'group = 1' ensures the line connects across factor levels within each facet
  geom_line(data = national_summary_2panel,
            aes(x = scenario, y = area_km2, group = 1, linetype = "National total"),
            linewidth = 1.2, color = "black") +
  geom_point(data = national_summary_2panel,
             aes(x = scenario, y = area_km2),
             size = 3, shape = 18, color = "black") +

  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +

  # Axis Scaling
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(
    name   = NULL,
    values = c("National total" = "dashed"),
    guide  = guide_legend(
      override.aes = list(color = "black", linewidth = 1.2, shape = NA)
    )
  ) +

  # Labels
  labs(
    x = NULL,
    y = expression(paste("Healthy coral area (", km^2, ")")),
    color = "Fishing Zone"
  ) +

  # Visual Theme
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(color = "grey40", size = 10, hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13)
  )

# Display the final plot
p_zone_area_2panel

ggsave(
  filename = file.path(out_dir_fig, "Fig_HealthyCoralArea_byFishing_zone.svg"),
  plot = p_zone_area_2panel,
  device = svglite::svglite,
  width = 12, height = 5.2, units = "in"
)

### Zone proportion line plot (two panels)
p_zone_prop_2panel <- ggplot(zones_healthy_2panel,
                              aes(x = scenario, y = percent_zone, group = zone, color = zone)) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = NULL,
    y = expression(paste("Proportion of Healthy coral area")),
    color = "Fishing Zone",
    caption = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(color = "grey40", size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13)
  )

p_zone_prop_2panel
ggsave(
  filename = file.path(out_dir_fig, "Fig_HealthyCoralProportion_byFishing_zone.svg"),
  plot = p_zone_prop_2panel,
  device = svglite::svglite,
  width = 12, height = 5.2, units = "in"
)


# Change detection map per zone
phase1_gain <- baseline_class == 1 & phase1_class == 2
phase2_gain <- baseline_class == 1 & phase2_class == 2
phase1_loss <- phase1_class == 2 & phase1_bleach_class == 1
phase2_loss <- phase2_class == 2 & phase2_bleach_class == 1

zone1_poly <- fishing_zones_utm[fishing_zones_utm$Name == "Zone 1", ]
gain_zone1 <- mask(
  phase2_improved,
  zone1_poly,
  overwrite = TRUE
)

plot(gain_zone1, main = "Phase 2 Coral Gains in Zone 1")

############################################
####
#### Describe pixels that recovery after phase 2
####
############################################

# upload tif files
# bathymetry
depth_res <- rast(
  file.path(out_dir, "depth_resampled_baseline.tif")
)

# baseline clasification
baseline_class <- rast(
  file.path(out_dir, "baseline_class.tif")
)

# phase2_class clasification
phase2_class <- rast(
  file.path(out_dir, "phase2_class.tif")
)

# baseline coral cover value
baseline <- rast(
  file.path(out_dir, "baseline_coral_cover.tif")
)

# phase 1 coral cover value
phase1 <- rast(
  file.path(out_dir, "phase1_coral_cover.tif")
)

phase1_bleach <- rast(
  file.path(out_dir, "bleach_coral_cover.tif")
)

bleaching_only <- rast(
  file.path(out_dir, "bleaching_only_coral_cover.tif")
)

phase2 <- rast(
  file.path(out_dir, "phase2_coral_cover.tif")
)


## Extract Coral Cover Values for Those Pixels
## extract values from each scenario raster only where recovered_pixels == TRUE
library(terra)
library(tidyverse)

# Extract depth
depth_vals <- values(depth_res)[recovered_cells]

#############
# Resample baseline to phase 1
template <- phase1  # grilla objetivo (10x10)

# 1) Resample baseline a la grilla de template
baseline_a <- resample(baseline, template, method = "bilinear")

# 2) Asegurar misma extensión: recortar todo a la intersección común
e_common <- intersect(ext(baseline_a), ext(phase1))
e_common <- intersect(e_common, ext(bleaching_only))
e_common <- intersect(e_common, ext(phase1_bleach))
e_common <- intersect(e_common, ext(phase2))

baseline_a        <- crop(baseline_a, e_common)
bleaching_only_a  <- crop(resample(bleaching_only, phase1, method="bilinear"), e_common)
phase1_a          <- crop(phase1, e_common)
phase1_bleach_a   <- crop(phase1_bleach, e_common)
phase2_a          <- crop(phase2, e_common)

# (Opcional) verifica que ahora sí calzan
compareGeom(baseline_a, phase1_a, stopOnError = TRUE)
compareGeom(phase1_a, phase1_bleach_a, stopOnError = TRUE)
compareGeom(phase1_a, phase2_a, stopOnError = TRUE)

# 3) Máscara de píxeles válidos
valid_mask <- !is.na(baseline_a) & !is.na(bleaching_only_a) & !is.na(phase1_a) & !is.na(phase1_bleach_a) & !is.na(phase2_a)
valid_cells <- which(values(valid_mask) == 1)

# 1) Píxeles válidos en TODOS los escenarios (evita NAs y longitudes distintas)
valid_mask <- !is.na(baseline) & !is.na(bleaching_only) & !is.na(phase1) & !is.na(phase1_bleach) & !is.na(phase2)

valid_cells <- which(values(valid_mask) == 1)

# ----------------------------
# 0) Align depth raster to Phase 1 grid + common extent
# ----------------------------
template <- phase1_a  # already cropped/aligned to Phase 1 grid in your workflow

depth_a <- resample(depth_res, template, method = "bilinear")

# crop depth to the exact extent (should already match, but keep it explicit)
depth_a <- crop(depth_a, ext(template))

stopifnot(compareGeom(depth_a, template, stopOnError = FALSE))

# ----------------------------
# 1) Valid pixels across ALL scenarios + depth
# ----------------------------
valid_mask <- !is.na(baseline_a) & !is.na(bleaching_only_a) & !is.na(phase1_a) & !is.na(phase1_bleach_a) & !is.na(phase2_a) & !is.na(depth_a)
valid_cells <- which(values(valid_mask) == 1)

# Extract vectors
baseline_v      <- values(baseline_a)[valid_cells]
bleaching_only_v <- values(bleaching_only_a)[valid_cells]
phase1_v        <- values(phase1_a)[valid_cells]
bleach_v        <- values(phase1_bleach_a)[valid_cells]
phase2_v        <- values(phase2_a)[valid_cells]

# Depth raster is negative (min -253, max 2). Usually depth is stored as negative.
# classify "shallow 0-10 m" using ABS(depth) in [0,10] and deep as >10.
depth_v <- values(depth_a)[valid_cells]

# ----------------------------
# 2) Build long table (all pixels) with depth group + dynamic cover class
# ----------------------------
trajectory_long <- tibble(
  pixel_id       = valid_cells,
  depth_m        = abs(depth_v),   # convert to positive meters (0 = surface)
  baseline       = baseline_v,
  bleaching_only = bleaching_only_v,
  phase1         = phase1_v,
  bleaching      = bleach_v,
  phase2         = phase2_v
) %>%
  mutate(
    depth_group = case_when(
      is.na(depth_m) ~ NA_character_,
      depth_m <= 10  ~ "Shallow (0–10 m)",
      depth_m > 10   ~ "Deep (>10 m)"
    ),
    depth_group = factor(depth_group, levels = c("Shallow (0–10 m)", "Deep (>10 m)"))
  ) %>%
  pivot_longer(
    cols = c(baseline, bleaching_only, phase1, bleaching, phase2),
    names_to = "scenario",
    values_to = "coral_cover"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = c("baseline", "bleaching_only", "phase1", "bleaching", "phase2")),
    coral_cover = pmin(pmax(coral_cover, 0), 100),
    cover_class = case_when(
      coral_cover < 5                       ~ "Critical (<5)",
      coral_cover >= 5  & coral_cover < 10  ~ "Poor (5–<10)",
      coral_cover >= 10 & coral_cover < 20  ~ "Fair (10–<20)",
      coral_cover >= 20 & coral_cover <= 40 ~ "Good (20–40)",
      coral_cover > 40                      ~ "Very good (>40)",
      TRUE                                  ~ NA_character_
    ),
    cover_class = factor(
      cover_class,
      levels = c("Critical (<5)", "Poor (5–<10)", "Fair (10–<20)", "Good (20–40)", "Very good (>40)")
    )
  ) %>%
  filter(!is.na(depth_group), !is.na(cover_class), !is.na(coral_cover))

# ----------------------------
# 3) Summarise mean + IQR by cover class, depth group, scenario
# ----------------------------
stats_depth_class <- trajectory_long %>%
  group_by(cover_class, depth_group, scenario) %>%
  summarise(
    n = n(),
    mean = mean(coral_cover, na.rm = TRUE),
    q25  = quantile(coral_cover, 0.25, na.rm = TRUE),
    q75  = quantile(coral_cover, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

stats_depth_class

# Plot (facets por cover class, líneas por profundidad, eje Y 0–50 + línea roja 15%)
depth_cols <- c(
  "Shallow (0–10 m)" = "#1F78B4",
  "Deep (>10 m)"     = "#33A02C"
)

ggplot(stats_depth_class,
       aes(x = scenario, y = mean, group = depth_group, color = depth_group, fill = depth_group)) +
  
  # 15% threshold line
  geom_hline(yintercept = 15, color = "#D62728", linewidth = 0.8) +
  
  # IQR band
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, color = NA) +
  
  # mean line + points
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  
  facet_wrap(~ cover_class, ncol = 3) +
  
  scale_x_discrete(labels = c("baseline" = "Baseline", "bleaching_only" = "Bleaching Only",
                              "phase1" = "Phase 1", "bleaching" = "Bleaching", "phase2" = "Phase 2")) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = depth_cols) +
  scale_fill_manual(values = depth_cols) +
  
  labs(
    title = "Coral cover trajectories by cover class and depth group",
    subtitle = "Lines = mean; shaded band = IQR (Q25–Q75). Facets are dynamic cover classes per scenario.",
    x = NULL,
    y = "Coral cover (%)",
    color = "Depth group",
    fill = "Depth group"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

#######################################################################################
## Direct test of the composition effect: Coral classes are defined statically 
## using baseline cover. By tracking these identical pixels across all scenarios, 
## group membership remains constant. This ensures that any observed changes 
## are interpreted without the bias of reclassification.
#######################################################################################

#######################################################################################
## Build long table using baseline_class_static + depth_group
#######################################################################################
library(dplyr)
library(tidyr)

# Build wide table for the SAME pixels across all scenarios
df_static <- tibble(
  pixel_id       = valid_cells,
  depth_m        = abs(values(depth_a)[valid_cells]),
  baseline       = values(baseline_a)[valid_cells],
  bleaching_only = values(bleaching_only_a)[valid_cells],
  phase1         = values(phase1_a)[valid_cells],
  bleaching      = values(phase1_bleach_a)[valid_cells],
  phase2         = values(phase2_a)[valid_cells]
) %>%
  mutate(
    depth_group = case_when(
      depth_m <= 10 ~ "Shallow (0–10 m)",
      depth_m > 10  ~ "Deep (>10 m)",
      TRUE          ~ NA_character_
    ),
    depth_group = factor(depth_group, levels = c("Shallow (0–10 m)", "Deep (>10 m)")),
    baseline_class_static = case_when(
      is.na(baseline)                 ~ NA_character_,
      baseline < 5                    ~ "Critical (<5)",
      baseline >= 5  & baseline < 10  ~ "Poor (5–<10)",
      baseline >= 10 & baseline < 20  ~ "Fair (10–<20)",
      baseline >= 20 & baseline <= 40 ~ "Good (20–40)",
      baseline > 40                   ~ "Very good (>40)"
    ),
    baseline_class_static = factor(
      baseline_class_static,
      levels = c("Critical (<5)", "Poor (5–<10)", "Fair (10–<20)", "Good (20–40)", "Very good (>40)")
    )
  ) %>%
  filter(!is.na(depth_group), !is.na(baseline_class_static))

trajectory_long_static <- df_static %>%
  pivot_longer(
    cols = c(baseline, bleaching_only, phase1, bleaching, phase2),
    names_to = "scenario",
    values_to = "coral_cover"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = c("baseline", "bleaching_only", "phase1", "bleaching", "phase2")),
    coral_cover = pmin(pmax(coral_cover, 0), 100)
  ) %>%
  filter(!is.na(coral_cover))

#######################################################################################
## Summary statistics by static class, depth, and scenario
#######################################################################################
stats_static <- trajectory_long_static %>%
  group_by(baseline_class_static, depth_group, scenario) %>%
  summarise(
    n = n(),
    mean = mean(coral_cover, na.rm = TRUE),
    q25  = quantile(coral_cover, 0.25, na.rm = TRUE),
    q75  = quantile(coral_cover, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

stats_static

# Plot por facets (cada panel = clase de baseline fija)
library(ggplot2)

depth_cols <- c(
  "Shallow (0–10 m)" = "#1F78B4",
  "Deep (>10 m)"     = "#33A02C"
)

p_facet <- ggplot(stats_static,
       aes(x = scenario, y = mean, group = depth_group, color = depth_group, fill = depth_group)) +
  geom_hline(yintercept = 15, color = "#D62728", linewidth = 0.8) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  
  facet_wrap(~ baseline_class_static, ncol = 3) +
  
  scale_x_discrete(labels = c("baseline" = "Baseline", "bleaching_only" = "Bleaching Only",
                              "phase1" = "Phase 1", "bleaching" = "Bleaching", "phase2" = "Phase 2")) +
  scale_y_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, 10),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_color_manual(values = depth_cols) +
  scale_fill_manual(values = depth_cols) +
  
  labs(
    title = NULL,
    subtitle = "Line = mean; band = IQR (Q25–Q75). Red line = 15% threshold.",
    x = NULL,
    y = "Coral cover (%)",
    color = "Depth group",
    fill = "Depth group"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    # Make facet separation obvious
    panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.7),
    panel.spacing = unit(0.9, "lines"),
    
    # Cleaner grids inside facets
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    # Facet strip styling (bigger + boxed)
    strip.background = element_rect(fill = "grey92", color = "grey25", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12),
    
    # Increase axis text size so it reads well within each facet
    axis.text.x = element_text(size = 11, color = "grey15"),
    axis.text.y = element_text(size = 11, color = "grey15"),
    axis.title.y = element_text(size = 12),
    
    # Titles
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

#save plot
# Guardar gráfico como PDF (Óptimo para Elsevier e Inkscape)
# Ajustamos los tamaños para que sean proporcionales al tamaño de exportación (180mm)
p_facet_export <- p_facet + 
  theme_minimal(base_size = 9) + # Bajamos la base de 13 a 9
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.subtitle = element_text(size = 8)
  )


ggsave(
  filename = file.path(out_dir, "Images_for_figure", "coral_cover_static_class_depth_facets_fixed.pdf"),
  plot = p_facet_export,
  device = cairo_pdf,
  width = 180 / 25.4,   # Ancho estándar Elsevier (Full page width)
  height = 100 / 25.4,  
  units = "in",
  dpi = 600
)


### two panels by depth
baseline_cols <- c(
  "Critical (<5)"   = "#D55E00",
  "Poor (5–<10)"    = "#E69F00",
  "Fair (10–<20)"   = "#009E73",
  "Good (20–40)"    = "#0072B2",
  "Very good (>40)" = "#CC79A7"
)

p_depth <- ggplot(
  stats_static,
  aes(
    x = scenario,
    y = mean,
    group = baseline_class_static,
    color = baseline_class_static,
    fill  = baseline_class_static
  )
) +
  
  geom_hline(yintercept = 15, color = "#D62728", linewidth = 0.8) +
  
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, color = NA) +
  
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  
  facet_wrap(~ depth_group, ncol = 2) +
  
  scale_x_discrete(labels = c("baseline" = "Baseline", "bleaching_only" = "Bleaching Only",
                              "phase1" = "Phase 1", "bleaching" = "Bleaching", "phase2" = "Phase 2")) +
  
  scale_y_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, 10),
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  scale_color_manual(values = baseline_cols) +
  scale_fill_manual(values = baseline_cols) +
  
  labs(
    title = NULL,
    subtitle = "Line = mean; band = IQR (Q25–Q75). Red line = 15% threshold.",
    x = NULL,
    y = "Coral cover (%)",
    color = "Baseline class",
    fill  = "Baseline class"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.7),
    panel.spacing = unit(0.9, "lines"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    strip.background = element_rect(fill = "grey92", color = "grey25", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12),
    
    axis.text.x = element_text(size = 11, color = "grey15"),
    axis.text.y = element_text(size = 11, color = "grey15"),
    axis.title.y = element_text(size = 12),
    
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

p_depth


############################################
# Coral recovery and fish biomass trajectories
# by baseline coral class and depth
#
# Output: 4-panel figure
# Rows  = depth (Shallow / Deep)
# Columns = metric (Coral cover / Fish biomass)
############################################


############################################
# 1. Load libraries
############################################

library(terra)
library(tidyverse)
library(patchwork)


############################################
# 2. Define directories
############################################

# Coral outputs
coral_dir <- file.path(out_dir, "spatial_files")

# Fish outputs
fish_dir <- "set/your/file/spatial_files"


############################################
# 3. Load coral rasters
############################################

depth_res <- rast(file.path(coral_dir, "depth_resampled_baseline.tif"))

baseline_coral       <- rast(file.path(coral_dir, "baseline_coral_cover.tif"))
bleaching_only_coral <- rast(file.path(coral_dir, "bleaching_only_coral_cover.tif"))
phase1_coral         <- rast(file.path(coral_dir, "phase1_coral_cover.tif"))
bleach_coral         <- rast(file.path(coral_dir, "bleach_coral_cover.tif"))
phase2_coral         <- rast(file.path(coral_dir, "phase2_coral_cover.tif"))


############################################
# 4. Load fish biomass rasters
############################################

baseline_fish <- rast(file.path(fish_dir, "baseline_fish_biomass_g100m2.tif"))
phase1_fish   <- rast(file.path(fish_dir, "phase1_fish_biomass_g100m2.tif"))
bleach_fish   <- rast(file.path(fish_dir, "phase1_bleach_fish_biomass_g100m2.tif"))
phase2_fish   <- rast(file.path(fish_dir, "phase2_fish_biomass_g100m2.tif"))


############################################
# 5. Align rasters to a common grid
############################################

template <- phase1_coral

baseline_coral_a       <- resample(baseline_coral,       template, method="bilinear")
bleaching_only_coral_a <- resample(bleaching_only_coral, template, method="bilinear")
baseline_fish_a        <- resample(baseline_fish,        template, method="bilinear")

e_common <- intersect(ext(baseline_coral_a), ext(phase1_coral))
e_common <- intersect(e_common, ext(bleaching_only_coral_a))
e_common <- intersect(e_common, ext(bleach_coral))
e_common <- intersect(e_common, ext(phase2_coral))

baseline_coral_a       <- crop(baseline_coral_a,       e_common)
bleaching_only_coral_a <- crop(bleaching_only_coral_a, e_common)
phase1_coral_a         <- crop(phase1_coral,           e_common)
bleach_coral_a         <- crop(bleach_coral,           e_common)
phase2_coral_a         <- crop(phase2_coral,           e_common)

baseline_fish_a  <- crop(baseline_fish_a, e_common)
phase1_fish_a    <- crop(phase1_fish, e_common)
bleach_fish_a    <- crop(bleach_fish, e_common)
phase2_fish_a    <- crop(phase2_fish, e_common)

depth_a <- resample(depth_res, template, method="bilinear")
depth_a <- crop(depth_a, e_common)


############################################
# 6. Identify valid pixels
############################################

valid_mask <-
  !is.na(baseline_coral_a) &
  !is.na(bleaching_only_coral_a) &
  !is.na(phase1_coral_a) &
  !is.na(bleach_coral_a) &
  !is.na(phase2_coral_a) &
  !is.na(baseline_fish_a) &
  !is.na(phase1_fish_a) &
  !is.na(bleach_fish_a) &
  !is.na(phase2_fish_a) &
  !is.na(depth_a)

valid_cells <- which(values(valid_mask) == 1)


############################################
# 7. Build base dataframe
############################################

df <- tibble(
  pixel_id = valid_cells,
  
  depth_m = abs(values(depth_a)[valid_cells]),
  
  coral_baseline       = values(baseline_coral_a)[valid_cells],
  coral_bleaching_only = values(bleaching_only_coral_a)[valid_cells],
  coral_phase1         = values(phase1_coral_a)[valid_cells],
  coral_bleach         = values(bleach_coral_a)[valid_cells],
  coral_phase2         = values(phase2_coral_a)[valid_cells],
  
  fish_baseline   = values(baseline_fish_a)[valid_cells],
  fish_phase1     = values(phase1_fish_a)[valid_cells],
  fish_bleach     = values(bleach_fish_a)[valid_cells],
  fish_phase2     = values(phase2_fish_a)[valid_cells]
)


############################################
# 8. Define depth groups
############################################

df <- df %>%
  mutate(
    depth_group = case_when(
      depth_m <= 10 ~ "Shallow (0–10 m)",
      depth_m > 10  ~ "Deep (>10 m)"
    )
  )


############################################
# 9. Define baseline coral classes
############################################

df <- df %>%
  mutate(
    baseline_class = case_when(
      coral_baseline < 5  ~ "Critical (<5)",
      coral_baseline < 10 ~ "Poor (5–<10)",
      coral_baseline < 20 ~ "Fair (10–<20)",
      coral_baseline <=40 ~ "Good (20–40)",
      coral_baseline >40  ~ "Very good (>40)"
    )
  )


############################################
# 10. Convert coral data to long format
############################################

coral_long <- df %>%
  select(depth_group, baseline_class,
         coral_baseline, coral_bleaching_only, coral_phase1, coral_bleach, coral_phase2) %>%
  pivot_longer(
    cols = starts_with("coral"),
    names_to="scenario",
    values_to="value"
  ) %>%
  mutate(
    metric="Coral cover"
  )


############################################
# 11. Convert fish data to long format
############################################

fish_long <- df %>%
  select(depth_group, baseline_class,
         fish_baseline, fish_phase1, fish_bleach, fish_phase2) %>%
  pivot_longer(
    cols = starts_with("fish"),
    names_to="scenario",
    values_to="value"
  ) %>%
  mutate(
    metric="Fish biomass"
  )


############################################
# 12. Merge datasets
############################################

data_long <- bind_rows(coral_long, fish_long)


############################################
# 13. Standardize scenario names
############################################

data_long <- data_long %>%
  mutate(
    scenario = case_when(
      str_detect(scenario,"bleaching_only") ~ "bleaching_only",
      str_detect(scenario,"baseline")       ~ "baseline",
      str_detect(scenario,"phase1")         ~ "phase1",
      str_detect(scenario,"bleach")         ~ "bleaching",
      str_detect(scenario,"phase2")         ~ "phase2"
    )
  )


############################################
# 14. Calculate statistics
############################################
library(patchwork)

stats <- stats %>%
  mutate(
    scenario = factor(scenario, 
                      levels = c("baseline", "bleaching_only", "phase1", "bleaching", "phase2")),
    baseline_class = factor(baseline_class,
                            levels = c("Critical (<5)", "Poor (5–<10)", "Fair (10–<20)", 
                                       "Good (20–40)", "Very good (>40)"))
  )

############################################
# 15. Color palette (color-blind safe)
############################################

baseline_cols <- c(
  "Critical (<5)"   = "#D55E00",
  "Poor (5–<10)"    = "#E69F00",
  "Fair (10–<20)"   = "#009E73",
  "Good (20–40)"    = "#0072B2",
  "Very good (>40)" = "#CC79A7"
)


############################################
# 16. Create final plot
############################################

# Plot CORAL COVER
p_coral <- ggplot(
  stats %>% filter(metric == "Coral cover"),
  aes(
    x = scenario,
    y = mean,
    group = baseline_class,
    color = baseline_class,
    fill = baseline_class
  )
) +
  
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.15,
    linewidth = 0.6,
    alpha = 0.7
  ) +
  
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  
  # coral threshold (SOLO EN CORAL)
  geom_hline(
    yintercept = 15,
    color = "#D62728",
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ depth_group,
    scales = "free_y",
    nrow = 1,
    ncol = 2
  ) +
  
  scale_x_discrete(
    labels = c(
      "baseline"       = "Baseline",
      "bleaching_only" = "Bleaching Only",
      "phase1"         = "Phase 1",
      "bleaching"      = "Bleaching",
      "phase2"         = "Phase 2"
    )
  ) +
  
  scale_color_manual(values = baseline_cols) +
  scale_fill_manual(values = baseline_cols) +
  
  labs(
    title = "Coral cover (%)",
    x = NULL,
    y = NULL,
    color = "Baseline coral class",
    fill = "Baseline coral class"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.7),
    panel.spacing = unit(0.9, "lines"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    strip.background = element_rect(fill = "grey92", color = "grey25", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12),
    
    axis.text.x = element_text(size = 11, color = "grey15"),
    axis.text.y = element_text(size = 11, color = "grey15"),
    
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    
    legend.position = "none"
  )


# plot without red line
p_fish <- ggplot(
  stats %>% filter(metric == "Fish biomass"),
  aes(
    x = scenario,
    y = mean,
    group = baseline_class,
    color = baseline_class,
    fill = baseline_class
  )
) +
  
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.15,
    linewidth = 0.6,
    alpha = 0.7
  ) +
  
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  
  # NO geom_hline
  
  facet_wrap(
    ~ depth_group,
    scales = "free_y",
    nrow = 1,
    ncol = 2
  ) +
  
  scale_x_discrete(
    labels = c(
      "baseline"       = "Baseline",
      "bleaching_only" = "Bleaching Only",
      "phase1"         = "Phase 1",
      "bleaching"      = "Bleaching",
      "phase2"         = "Phase 2"
    )
  ) +
  
  # same scale for both biomass panels 0-4000 
  scale_y_continuous(limits = c(0, 4000)) +
  
  scale_color_manual(values = baseline_cols) +
  scale_fill_manual(values = baseline_cols) +
  
  labs(
    title = "Fish biomass (g 100m⁻²)",
    x = NULL,
    y = NULL,
    color = "Baseline coral class",
    fill = "Baseline coral class"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.7),
    panel.spacing = unit(0.9, "lines"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    strip.background = element_rect(fill = "grey92", color = "grey25", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12),
    
    axis.text.x = element_text(size = 11, color = "grey15"),
    axis.text.y = element_text(size = 11, color = "grey15"),
    
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )


############################################
# 17. Display plot
############################################
p <- p_coral / p_fish

p
