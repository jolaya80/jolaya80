

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
    cols = c(baseline, phase1, bleaching, phase2),
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

# Opcional: Reordenar los escenarios para que tengan sentido lógico
fishery_estimates$scenario <- factor(fishery_estimates$scenario, 
                                 levels = c("baseline", "phase1", "bleaching", "phase2"))

ggplot(fishery_estimates, aes(x = scenario, y = estimated_catch_tons, fill = Location)) +
  # Dibujar barras con posición lado a lado
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  # Añadir etiquetas de texto sobre las barras
  geom_text(
    aes(label = round(estimated_catch_tons, 2)), 
    position = position_dodge(width = 0.9), 
    vjust = -0.5,     # Ajuste vertical para que esté arriba de la barra
    size = 3.5,       # Tamaño de la fuente
    fontface = "bold"
  ) +
  # Paleta de colores profesional
  scale_fill_brewer(palette = "Set1") +
  # Títulos y etiquetas en inglés
  labs(
    x = "Scenario",
    y = "Annual Catch (tons)",
    fill = "Location"
  ) +
  # Ajustar límites del eje Y para que el texto no se corte
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  # Estética limpia
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(), # Quitar líneas verticales para más limpieza
    plot.title = element_text(face = "bold", size = 14)
  )

## line plot
ggplot(fishery_estimates, aes(x = scenario, y = estimated_catch_tons, 
                              color = Location, group = Location)) +
  # Add the lines (using a slightly thicker line for better visibility)
  geom_line(linewidth = 1.2) + 
  # Add points to clearly mark each scenario
  geom_point(size = 3) +
  # Add text labels on top of the points
  geom_text(
    aes(label = round(estimated_catch_tons, 2)), 
    vjust = -1.2,     # Adjust vertically to sit above the points
    size = 3.5,       
    fontface = "bold",
    show.legend = FALSE # Prevents 'a' appearing in the legend
  ) +
  # Professional color palette for lines (Color instead of Fill)
  scale_color_brewer(palette = "Set1") +
  # Titles and labels
  labs(
    x = "Scenario",
    y = "Annual Catch (tons)",
    fill = "Location"
  ) +
  # Adjust Y-axis to prevent labels from being cut off
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  # Clean aesthetics
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold")
  )



















