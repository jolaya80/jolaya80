
############################################################
### CPUE estimation & spatial pipeline
############################################################

library(tidyverse)
library(performance)
library(car)
library(rsq)
library(ggplot2)
library(ggpmisc)
library(sf)
library(raster)
library(terra)
library(tools)
library(forcats)

## =========================
## 1. Load data
## =========================

data <- read_csv(
  "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/data/Fishery_surveys/fish_catch_bz.csv",
  show_col_types = FALSE
)

data <- data %>%
  mutate(
    # BUG ANNOTATION: Verify column naming - assumes catch columns all start with "catch_"
    total_catch_lb = rowSums(across(starts_with("catch_"), .names = "catch_"), na.rm = TRUE),
    cpue_raw_hours = total_catch_lb / hours_trip_fish  # Standardized by fishing hours
  )

# Zero values treated as missing (indicate non-response, not actual zero catch)
data <- data %>%
  mutate(
    cpue_raw_hours = na_if(cpue_raw_hours, 0)
  )

## =========================
## 2. Clean categorical variables (Location, fishing_type_main)
## =========================
# Extract primary categories from ranked responses and set factor levels

# extract first category from the ranking questions
data <- data %>%
  mutate(
    fishing_type_main = str_split(fishing_type_rank, ",") %>% map_chr(1),
    main_buyer        = str_split(main_buyers_fish, ",")  %>% map_chr(1)
  )

# define levels
lvl_loc     <- data$Location          %>% na.omit() %>% unique() %>% sort()
lvl_fishing <- data$fishing_type_main %>% na.omit() %>% unique() %>% sort()
lvl_buyer   <- data$main_buyer        %>% na.omit() %>% unique() %>% sort()
lvl_sea_adapt <- data$Seasonal_Adaptability %>% na.omit() %>% unique %>% sort()

# convert to factors with fixed levels
data <- data %>%
  mutate(
    Location          = factor(Location,          levels = lvl_loc),
    fishing_type_main = factor(fishing_type_main, levels = lvl_fishing),
    main_buyer        = factor(main_buyer, levels = lvl_buyer),
    Seasonal_Adaptability = factor(Seasonal_Adaptability, levels = lvl_sea_adapt)
  )

summary(data)

## =========================
## 3. Exploratory data visualization
## =========================

# Relationship between fishing hours and raw CPUE
ggplot(data, aes(x = hours_trip_fish, y = cpue_raw_hours)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Raw CPUE vs Hours per Trip",
       x = "Hours per Trip",
       y = "CPUE (lb/hour)") +
  theme_bw()

# Create professional boxplot with individual data points showing CPUE distribution
plot_cpue <- ggplot(data, aes(x = fishing_type_main, y = cpue_raw_hours, fill = fishing_type_main)) +
  
  # Boxplot layer showing quartiles and median
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  
  # Individual data points with jitter for visibility
  geom_jitter(
    aes(color = fishing_type_main),
    width = 0.15,
    size = 1,
    alpha = 0.8
  ) +
  
  # Separate panels by Location
  facet_wrap(~Location) +
  
  # Color scheme
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  
  # Labels and title
  labs(
    x = NULL,
    y = "CPUE (lb/hour)",
    title = NULL
  ) +
  
  # Theme settings
  theme_classic(base_size = 10, base_family = "serif") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    # X-axis angle adjustment: 45 degrees with right alignment
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_cpue

# save as PDF (Vectorial)
#define output folder
target_folder <- "G:/Shared drives/NSF CoPE internal/2 - Deliverables/Publications/Olaya_et_al_Belize_FisheryModel/figures"

# Create full file path combining folder and filename
file_path <- file.path(target_folder, "supp_Figure_CPUE_fishing_type.tif")

# Save with journal publication specifications
ggsave(
  filename = file_path,
  plot = plot_cpue,
  device = "tiff",
  width = 85,           # Single-column width (mm)
  height = 90,          # Proportional height
  units = "mm",
  dpi = 600,            # High resolution for line and point combinations
  compression = "lzw"   # Lossless compression required by most journals
)

## =========================
## 4. GLM-based CPUE standardization and model selection
## =========================
# RESEARCH OBJECTIVE: Standardize CPUE across fisher populations to account for 
# differences in fishing practices, location, and socioeconomic characteristics.
# This allows for fair comparison of fishing productivity across space and time.

# WHY GAMMA GLM WITH LOG-LINK?
# CPUE data are positive, continuous, and typically right-skewed (many low values, 
# few high values). The Gamma distribution natural handles this asymmetry better 
# than normal distribution. The log-link ensures predicted values remain positive.

# VARIABLES AND HYPOTHESES:
# - Location: Different zones may have different fish availability/accessibility
# - Seasonal_Adaptability: Specialists vs generalists may differ in efficiency
# - fishing_type_main: Different gear/methods affect catch rates
# - coop_membership: Group membership may affect market access or sharing practices
# - main_buyer: Buyer type may affect effort allocation or catch composition

# Variables required across all candidate models
vars_full <- c(
  "cpue_raw_hours",
  "Location",
  "Seasonal_Adaptability",
  "fishing_type_main",
  "coop_membership",
  "main_buyer"
)

# Complete cases only (ensures valid model comparison)
# Removes rows with missing values in any predictor or outcome
data_cc <- data[complete.cases(data[, vars_full]), ] %>%
  mutate(
    Location = factor(Location),
    fishing_type_main = factor(fishing_type_main),
    Seasonal_Adaptability = factor(Seasonal_Adaptability),
    coop_membership = factor(coop_membership, levels = c("No", "Yes")),
    main_buyer = factor(main_buyer)
  )

# MODEL 1: Full model with all predictors
# Formula: log(CPUE) = Location + Seasonal_Adaptability + fishing_type + coop_membership + buyer_type
# Captures all measured socioeconomic and ecological factors
m_full_cc <- glm(
  cpue_raw_hours ~ Location + Seasonal_Adaptability +
    fishing_type_main + coop_membership + main_buyer,
  family = Gamma(link = "log"),
  data = data_cc
)

# MODEL 2: Reduced model - main ecological/behavioral effects
# Removes buyer and coop membership (socioeconomic factors)
# Tests whether simpler model adequately captures CPUE variation
m_red1_cc <- glm(
  cpue_raw_hours ~ Location + Seasonal_Adaptability +
    fishing_type_main,
  family = Gamma(link = "log"),
  data = data_cc
)

# MODEL 3: Minimal model - only location and fishing type
# Maximum parsimony; primarily for robustness checks
# Tests whether spatial and gear effects alone explain variation
m_red2_cc <- glm(
  cpue_raw_hours ~ Location + fishing_type_main,
  family = Gamma(link = "log"),
  data = data_cc
)

# Detailed summary table
summary_table <- broom::tidy(m_red2_cc) %>%
  mutate(
    exp_beta = exp(estimate)  # This calculates the "fold-change" / Exp(beta)
  )

print(summary_table)

# MODEL COMPARISON USING INFORMATION CRITERIA AND HYPOTHESIS TESTS
# AIC (Akaike Information Criterion) balances model fit with parsimony
# Lower AIC indicates better predictive performance relative to model complexity
# Likelihood ratio tests compare nested models:
#   p > 0.05 suggests simpler model is adequate (dropped variables not significant)
#   p < 0.05 suggests more complex model needed (dropped variables matter)

AIC(m_full_cc, m_red1_cc, m_red2_cc)

# Test whether coop_membership and main_buyer improve model fit
anova(m_full_cc, m_red1_cc, test = "Chisq")

# Test whether Seasonal_Adaptability improves model fit
anova(m_red1_cc, m_red2_cc, test = "Chisq")

# Test whether Seasonal_Adaptability improves model fit
anova(m_full_cc, m_red2_cc, test = "Chisq")

## =========================
## 5. Model diagnostics
## =========================

coef_df <- bind_rows(
  broom::tidy(m_full_cc, conf.int = TRUE)  %>% mutate(Model = "Full (all predictors)"),
  broom::tidy(m_red1_cc, conf.int = TRUE) %>% mutate(Model = "Reduced"),
  broom::tidy(m_red2_cc, conf.int = TRUE) %>% mutate(Model = "Minimun")
) %>%
  filter(term != "(Intercept)")

ggplot(coef_df,
       aes(x = estimate, y = term, color = Model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Log-scale effect on CPUE",
    y = NULL,
    title = "Comparison of coefficient estimates across model specifications"
  ) +
  theme_bw()

## =========================
## 8. Predicted CPUE for a reference fisher (intuitive comparison)
## =========================

# reference (baseline) fisher
newdata <- expand.grid(
  Location = levels(data_cc$Location),
  fishing_type_main = levels(data_cc$fishing_type_main),
  Seasonal_Adaptability = "Fish_specialist",
  coop_membership = "No",
  main_buyer = "Community"
)

# predictions (only variables used by each model)
pred_df <- bind_rows(
  
  newdata %>%
    mutate(CPUE = predict(m_full_cc, newdata, type = "response"),
           Model = "Full"),
  
  newdata %>%
    mutate(CPUE = predict(m_red1_cc, newdata, type = "response"),
           Model = "Reduced"),
  
  newdata %>%
    mutate(CPUE = predict(m_red2_cc, newdata, type = "response"),
           Model = "Minimal")
)

ggplot(pred_df,
       aes(x = fishing_type_main, y = CPUE, fill = Location)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Model) +
  labs(
    y = "Predicted CPUE (lb per hour)",
    x = "Fishing strategy"
  ) +
  theme_bw()

## Model validation: Compare predicted vs observed CPUE

# Generate predictions for all observed fishers
data_pred <- data_cc %>%
  mutate(
    CPUE_pred = predict(m_red2_cc, newdata = data_cc, type = "response")
  )

# Calculate correlation between predicted and observed
cor_val <- cor(data_pred$CPUE_pred, data_pred$cpue_raw_hours, use = "complete.obs")
cor_val

# Scatterplot of observed vs predicted CPUE
# Points along the 1:1 line (red dashed) indicate perfect prediction
ggplot(data_pred, aes(x = cpue_raw_hours , y = CPUE_pred)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  labs(
    x = "Observed CPUE (lb per hour)",
    y = "Predicted CPUE (lb per hour)",
    title = "Model validation: Observed vs Predicted CPUE"
  ) +
  theme_bw(base_size = 12) + 
  annotate("text", x = max(data_pred$CPUE_pred, na.rm=TRUE)*0.8, 
           y = max(data_pred$cpue_raw_hours, na.rm=TRUE)*0.9,
           label = paste0("Correlation (r) = ", round(cor_val, 2)), size = 4, fontface="bold")

## =========================
## 9. Standardization of CPUE (hours-based model)
## =========================

# Creamos el boxplot profesional con puntos
ggplot(data_cc, aes(x = Location, y = cpue_raw_hours, fill = fishing_type_main)) +
  
  # Boxplot
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  
  # Puntos individuales, NA (si hubiera) como círculos vacíos
  geom_jitter(aes(color = fishing_type_main, shape = is.na(cpue_raw_hours)),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, alpha = 0.8) +
  
  # Formas y colores
  scale_shape_manual(values = c(16, 1), labels = c("Value", "NA")) +  # 16 = punto lleno, 1 = círculo vacío
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  
  # Etiquetas
  labs(
    x = "Location",
    y = "CPUE (lb/h)",
    fill = "fishing_type_main",
    color = "fishing_type_main",
    shape = "Data Type",
    title = "Distribution of CPUE by Location and fishing_type_main"
  ) +
  
  # Tema limpio
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


## =========================
## 10. Export CPUE table
## =========================

Table_output_folder <- "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results"

if (!dir.exists(Table_output_folder)) {
  dir.create(Table_output_folder, recursive = TRUE)
}

write.csv(
  data_pred,
  file      = file.path(Table_output_folder, "finfish_survey_data_cpue.csv"),
  row.names = FALSE
)

## =========================
## 12. Join CPUE table to fisher polygons
## =========================

# read data with effort
effort <- read_csv(
  "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/cpue_data_clean.csv",
  show_col_types = FALSE
)

# Unimos solo la columna deseada usando globalid como llave
data_pred <- data_pred %>%
  left_join(effort %>% dplyr::select(globalid, months_fishing), by = "globalid")

# Verificamos que se haya añadido al final
summary(data_pred$months_fishing)

cpue_data <- data_pred %>%
  dplyr::select(
    globalid,
    Location,
    fishingzone_licensed,
    Seasonal_Adaptability,
    gear_finfish,
    months_fishing,         # Mᵢ
    number_trips_fish,      # Tᵢ
    hours_trip_fish,        # Hᵢ
    total_catch_lb,         # Cᵢ
    CPUE_pred,
    fishing_type_main,
    number_trips_fish,
    hours_trip_fish
  )


fish_polygons <- st_read(
  "G:/Shared drives/NSF CoPE internal/GIS_CoPE/GIS_Belize/2_model_inputs_belize/Fishery_interviews/shp/fish_polygons_fishers.shp"
)

fishing_zones_utm <- st_read(
  "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/01_shp/fishing_zones_utm.gpkg"
)

#######################################################################################################################

#########################################################################################

# Clean data: Keep fishers with valid CPUE estimates
# Note: Missing zone information is retained as NA
data_pred <- data_pred %>%
  mutate(
    annual_effort = if_else(!is.na(hours_trip_fish) & !is.na(months_fishing) & !is.na(number_trips_fish),
                            months_fishing * number_trips_fish * hours_trip_fish, NA_real_),
    annual_catch = if_else(!is.na(total_catch_lb) & !is.na(months_fishing) & !is.na(number_trips_fish),
                           months_fishing * number_trips_fish * total_catch_lb, NA_real_)
  ) %>%
  filter(!is.na(CPUE_pred) | !is.na(annual_effort))

write.csv(
  data_pred,
  "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/cpue_predicted_GLM.csv",
  row.names = FALSE
)

# CPUE Location cpue_location <- cpue_data_clean %>%
cpue_location <- cpue_data %>%
  group_by(Location) %>%
  summarise(
    total_annual_catch  = sum(annual_catch, na.rm = TRUE),
    total_annual_effort = sum(annual_effort, na.rm = TRUE),
    cpue_annual         = total_annual_catch / total_annual_effort,
    n_fishers           = n(),
    .groups = "drop"
  )


# cpue by Location x type fisher
cpue_location_type <- cpue_data_clean %>%
  group_by(Location, fishing_type_main) %>%
  summarise(
    total_annual_catch  = sum(annual_catch, na.rm = TRUE),
    total_annual_effort = sum(annual_effort, na.rm = TRUE),
    cpue_annual         = total_annual_catch / total_annual_effort,
    n_fishers = n(),
    .groups = "drop"
  )

write.csv(
  cpue_location_type,
  "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/cpue_location_type.csv",
  row.names = FALSE
)



# bring in biomass per zone (all scenarios)
biomass_zone <- read_csv(
  "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/Fish_biomass/outputs/total_biomass_zone_scenarios.csv"
)

# compute scenario-specific CPUE per fisher
cpue_fisher_scenarios <- cpue_fisher_zoned %>%
  filter(
    !is.na(cpue_raw_hours),
    !is.na(annual_effort)
  ) %>%
  left_join(
    biomass_zone,
    by = "zone_id"
  ) %>%
  group_by(globalid) %>%
  mutate(
    # Extract baseline biomass ONLY if zone exists
    biomass_baseline = if_else(
      !is.na(zone_id),
      Biomass_ton[Scenario == "Baseline"],
      NA_real_
    ),
    
    # Scale CPUE ONLY when biomass is known
    cpue_scenario = if_else(
      !is.na(zone_id),
      cpue_raw_hours * (Biomass_ton / biomass_baseline),
      NA_real_
    ),
    
    # Compute scenario catch ONLY when scaling is possible
    catch_scenario = if_else(
      !is.na(zone_id),
      cpue_scenario * annual_effort,
      NA_real_
    )
  ) %>%
  ungroup()

cpue_fisher_scenarios <- cpue_fisher_scenarios %>%
  mutate(
    needs_manual_zone = is.na(zone_id)
  )

write.csv(cpue_fisher_scenarios,
          "C:/Users/jolaya/Documents/GitHub_projects/Networks_SSF_NatCap/models/CPUE_fisheries/results/cpue_fisher_scenarios.csv",
          row.names = FALSE)

fishers_no_zone <- cpue_fisher_scenarios %>%
  filter(needs_manual_zone) %>%
  dplyr::select(globalid, fishingzone_licensed, Location, fishing_type_main, cpue_raw_hours)



# NOTE: Direct zone assignment to fishers SKIPPED as per user request

########## Aggregate fisher results to zone × scenario
zone_catch_scenarios <- cpue_fisher_scenarios %>%
  group_by(zone_id_final, Zone.y, Scenario) %>%
  summarise(
    total_catch = sum(catch_scenario, na.rm = TRUE),
    total_effort = sum(annual_effort, na.rm = TRUE),
    cpue_zone = total_catch / total_effort,
    .groups = "drop"
  )

# Compute biomass & catch change relative to baseline
zone_baseline <- zone_catch_scenarios %>%
  filter(Scenario == "Baseline") %>%
  dplyr::select(
    zone_id_final,
    cpue_baseline = cpue_zone,
    catch_baseline = total_catch
  )

# Join and compute proportional change
zone_change <- zone_catch_scenarios %>%
  left_join(zone_baseline, by = "zone_id_final") %>%
  left_join(
    biomass_zone %>%
      filter(Scenario == "Baseline") %>%
      dplyr::select(zone_id, biomass_baseline = Biomass_ton),
    by = c("zone_id_final" = "zone_id")
  ) %>%
  left_join(
    biomass_zone,
    by = c("zone_id_final" = "zone_id", "Scenario")
  ) %>%
  mutate(
    biomass_change = Biomass_ton / biomass_baseline,
    catch_change   = total_catch / catch_baseline
  )

# plot library(ggplot2)

ggplot(
  zone_change %>% filter(Scenario != "Baseline"),
  aes(x = biomass_change, y = catch_change, color = Scenario)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 4) +
  geom_text(
    aes(label = Zone),
    vjust = -0.8,
    size = 3.5
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Relative Fish Biomass (Scenario / Baseline)",
    y = "Relative Catch (Scenario / Baseline)",
    title = "Catch response to biomass change across fishing zones",
    subtitle = "Effort-weighted aggregation of fisher CPUE under restoration and bleaching scenarios",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)

## other plot
zone_biomass_catch <- zone_catch_scenarios %>%
  left_join(
    biomass_zone,
    by = c("zone_id_final" = "zone_id", "Scenario")
  ) %>%
  filter(!is.na(Biomass_ton), !is.na(total_catch))

library(ggplot2)

ggplot(
  zone_biomass_catch,
  aes(
    x = Biomass_ton,
    y = total_catch,
    color = Scenario
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = Zone),
    vjust = -0.7,
    size = 3.5
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Total Fish Biomass (tons)",
    y = "Total Annual Catch (lb/year)",
    title = "Fish Catch as a Function of Biomass Across Fishing Zones",
    subtitle = "Zone-level aggregation of fisher catches under restoration and bleaching scenarios",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)






