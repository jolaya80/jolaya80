

############################################################
### SCENARIO ANALYSIS – Expected Catch under Biomass Change
### Effort constant, CPUE scales with biomass
############################################################

library(tidyverse)

#-----------------------------------------------------------
# 1. Read data
#-----------------------------------------------------------

biomass_path <- "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/biomass_by_scenarios_tons.csv"

cpue_path <- "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/cpue_predicted_GLM.csv"

catch_path <- "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/finfish_survey_data_cpue.csv"

fish_biomass_raw <- read_csv(biomass_path)
cpue_data        <- read_csv(cpue_path)
catch_data       <- read_csv(catch_path)

#-----------------------------------------------------------
# 2. Keep only Zone 1 (Caye Caulker) and Zone 3 (Placencia)
#-----------------------------------------------------------

zone_lookup <- tibble(
  zone_id  = c(1, 3),
  Location = c("Caye Caulker", "Placencia")
)

fish_biomass <- fish_biomass_raw %>%
  filter(zone_id %in% c(1, 3)) %>%
  left_join(zone_lookup, by = "zone_id")

#-----------------------------------------------------------
# 3. Convert biomass to long format and compute ratios
#-----------------------------------------------------------

biomass_long <- fish_biomass %>%
  pivot_longer(
    cols = c(baseline, phase1, bleaching_only, bleaching_only_phase2, bleaching, phase2),
    names_to = "scenario",
    values_to = "biomass_tons"
  )

# Extract baseline biomass per zone
baseline_vals <- biomass_long %>%
  filter(scenario == "baseline") %>%
  select(zone_id, baseline_biomass = biomass_tons)

biomass_long <- biomass_long %>%
  left_join(baseline_vals, by = "zone_id") %>%
  mutate(
    biomass_ratio = biomass_tons / baseline_biomass
  )

#-----------------------------------------------------------
# 4. Use ONLY baseline CPUE per locality (from GLM output)
#-----------------------------------------------------------

# Calculate CPUE based on location-specific baselines
biomass_with_cpue <- biomass_long %>%
  mutate(
    # Assign the correct baseline CPUE based on the Location
    baseline_cpue = case_when(
      Location == "Caye Caulker" ~ 5.57,
      Location == "Placencia"    ~ 39.3,
      TRUE ~ NA_real_
    ),
    # Calculate the scenario CPUE using the proportional change (biomass_ratio)
    cpue_lb_h = baseline_cpue * biomass_ratio
  )

# Review the results
biomass_with_cpue %>%
  dplyr::select(Location, scenario, biomass_ratio, baseline_cpue, cpue_lb_h)

#-----------------------------------------------------------
# 5. Compute mean annual effort per locality
#-----------------------------------------------------------

effort_summary <- cpue_data %>%
  filter(Location %in% c("Caye Caulker", "Placencia")) %>%
  group_by(Location) %>%
  summarise(
    mean_annual_effort = mean(annual_effort, na.rm = TRUE),
    .groups = "drop"
  )

#-----------------------------------------------------------
# 6. MAIN ESTIMATION
#-----------------------------------------------------------
# Logic:
# CPUE_scenario = CPUE_baseline × biomass_ratio
# Catch(lb/year) = CPUE(lb/h) × Effort(h/year)

LBS_PER_TON <- 2204.62

# Join tables and calculate metrics
fishery_estimates <- biomass_with_cpue %>%
  left_join(effort_summary, by = "Location") %>%
  mutate(
    # Percentage change in CPUE relative to baseline
    cpue_pct_change = (biomass_ratio - 1) * 100,
    
    # Estimated total catch (CPUE * Effort)
    estimated_catch_lb = cpue_lb_h * mean_annual_effort,
    
    # Convert total catch from pounds to metric tons
    estimated_catch_tons = estimated_catch_lb / LBS_PER_TON
  )

#-----------------------------------------------------------
# 7. Final results table
#-----------------------------------------------------------

# 3. View the final results
fishery_estimates %>%
  dplyr::select(Location, scenario, cpue_lb_h, cpue_pct_change, mean_annual_effort, estimated_catch_lb, estimated_catch_tons)

print(fishery_estimates)

write.csv(fishery_estimates, "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/results_table.csv", row.names = FALSE)

library(ggplot2)
library(scales)
library(forcats)

# Reorder scenarios in logical temporal order (including new No Intervention scenarios)
fishery_estimates$scenario <- factor(fishery_estimates$scenario,
                                     levels = c("baseline", "phase1",
                                                "bleaching_only", "bleaching_only_phase2",
                                                "bleaching", "phase2"),
                                     labels = c("Baseline", "Phase 1",
                                                "Bleaching Only", "Bleaching Only Phase 2",
                                                "Bleaching", "Phase 2"))

#-----------------------------------------------------------
# Two-panel comparative structures
# Mirrors the pattern in fish_biomass_scenarios_calculation.R
#-----------------------------------------------------------

## 3c. Two-panel data: No Intervention vs With Restoration
fishery_estimates_2panel <- bind_rows(
  fishery_estimates %>%
    filter(scenario %in% c("Baseline", "Bleaching Only")) %>%
    mutate(scenario_group = "No Intervention"),
  fishery_estimates %>%
    filter(scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(scenario_group = "With Restoration")
) %>%
  mutate(scenario_group = factor(scenario_group,
                                 levels = c("No Intervention", "With Restoration")))

## 3d. Two-panel data with temporal x-axis labels
fishery_estimates_2panel_temporal <- bind_rows(

  # ---- No Intervention panel (4 synthetic time points) ----
  fishery_estimates %>%
    filter(scenario == "Baseline") %>%
    mutate(scenario = "Year zero", scenario_group = "No Intervention"),
  fishery_estimates %>%
    filter(scenario == "Baseline") %>%
    mutate(scenario = "1 Year", scenario_group = "No Intervention"),
  fishery_estimates %>%
    filter(scenario == "Bleaching Only") %>%
    mutate(scenario = "2 Years", scenario_group = "No Intervention"),
  fishery_estimates %>%
    filter(scenario == "Bleaching Only Phase 2") %>%
    mutate(scenario = "~5 Years", scenario_group = "No Intervention"),

  # ---- With Restoration panel (renamed to temporal labels) ----
  fishery_estimates %>%
    filter(scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(
      scenario = recode(as.character(scenario),
                        "Baseline"  = "Year zero",
                        "Phase 1"   = "1 Year",
                        "Bleaching" = "2 Years",
                        "Phase 2"   = "~5 Years"),
      scenario_group = "With Restoration"
    )

) %>%
  mutate(
    scenario       = factor(scenario, levels = c("Year zero", "1 Year", "2 Years", "~5 Years")),
    scenario_group = factor(scenario_group, levels = c("No Intervention", "With Restoration")),
    Location       = factor(Location)
  )

#-----------------------------------------------------------
# Location colour palette (matches Dark2 style used in biomass script)
#-----------------------------------------------------------
loc_levels <- levels(fishery_estimates_2panel_temporal$Location)
loc_cols   <- setNames(
  RColorBrewer::brewer.pal(n = max(3, length(loc_levels)), name = "Dark2")[seq_along(loc_levels)],
  loc_levels
)

#-----------------------------------------------------------
# 3e. Two-panel bar plot: No Intervention vs With Restoration
#-----------------------------------------------------------
p_catch_bar_2panel <- ggplot(fishery_estimates_2panel,
                              aes(x = scenario, y = estimated_catch_tons, fill = Location)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = round(estimated_catch_tons, 2)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    x     = NULL,
    y     = "Annual Catch (tons)",
    fill  = "Location"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(
    legend.position      = "top",
    panel.grid.major.x   = element_blank(),
    plot.title           = element_text(face = "bold", size = 14),
    strip.text           = element_text(face = "bold", size = 13),
    axis.text.x          = element_text(angle = 20, hjust = 1)
  )

p_catch_bar_2panel

#-----------------------------------------------------------
# National total per scenario_group (dashed overlay line for line plot)
#-----------------------------------------------------------
national_catch_2panel <- fishery_estimates_2panel_temporal %>%
  group_by(scenario_group, scenario) %>%
  summarise(estimated_catch_tons = sum(estimated_catch_tons, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Location = "National")

#-----------------------------------------------------------
# 3f. Two-panel line plot with temporal x-axis labels
#-----------------------------------------------------------
p_catch_line_2panel <- ggplot(fishery_estimates_2panel_temporal,
                               aes(x = scenario, y = estimated_catch_tons,
                                   group = Location, color = Location)) +

  # Regional lines
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +

  # National total (dashed black line)
  geom_line(data = national_catch_2panel,
            aes(x = scenario, y = estimated_catch_tons, group = 1,
                linetype = "National total"),
            linewidth = 1.2, color = "black") +
  geom_point(data = national_catch_2panel,
             aes(x = scenario, y = estimated_catch_tons),
             size = 3, shape = 18, color = "black") +

  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +

  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  scale_color_manual(values = loc_cols) +
  scale_linetype_manual(
    name   = NULL,
    values = c("National total" = "dashed"),
    guide  = guide_legend(
      override.aes = list(color = "black", linewidth = 1.2, shape = NA)
    )
  ) +
  labs(
    x     = NULL,
    y     = "Annual Catch (tons)",
    color = "Location"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position  = "plot",
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    axis.text.x          = element_text(angle = 20, hjust = 1),
    axis.title           = element_text(face = "bold"),
    legend.position      = "right",
    strip.text           = element_text(face = "bold", size = 13)
  )

p_catch_line_2panel

#-----------------------------------------------------------
# 3g. Save updated outputs
#-----------------------------------------------------------
results_folder <- "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results"

write.csv(fishery_estimates,
          file.path(results_folder, "results_table.csv"),
          row.names = FALSE)

ggsave(
  filename = file.path(results_folder, "Fig_catch_tons_2panel_bar.svg"),
  plot     = p_catch_bar_2panel,
  device   = svglite::svglite,
  width    = 12,
  height   = 5.2,
  units    = "in",
  fix_text_size = FALSE
)

ggsave(
  filename = file.path(results_folder, "Fig_catch_tons_2panel_line.svg"),
  plot     = p_catch_line_2panel,
  device   = svglite::svglite,
  width    = 12,
  height   = 5.2,
  units    = "in",
  fix_text_size = FALSE
)



















