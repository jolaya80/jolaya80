# ============================================================
# Script 01: Fisher Bioeconomic Model
# Project:   Coral Restoration Benefits — Belize SSF
# Purpose:   Multispecies catch partitioning → social allocation
#            → economic valuation → food security metrics
#            → scenario delta calculations
# Author:    [Your name]
# Date:      2026
# ============================================================

setwd("~/Documents/GitHub_projects/Networks_SSF_NatCap/models/Economic")

# ============================================================
# 0. Packages
# ============================================================
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(forcats)

# ============================================================
# 1. Constants
# ============================================================
TON_TO_LB   <- 2204.62   # 1 metric ton = 2204.62 lb
LB_PER_MEAL <- 0.5       # 1 meal ≈ 0.5 lb (227 g) per adult portion
                          # (FAO/fisheries nutrition literature standard)

# Scenario display order (all 6 scenarios across both panels)
SCENARIO_ORDER <- c("Baseline", "Phase 1", "Bleaching", "Phase 2",
                    "Bleaching Only", "Bleaching Only Phase 2")

# ============================================================
# 2. Load data
# ============================================================

# 2a. Catch & biomass by scenario (PRIMARY CATCH INPUT)
# Reads from catch_estimation.R output (results_table.csv) which includes all
# 6 scenarios: Baseline, Phase 1, Bleaching, Phase 2 (With Restoration) plus
# Bleaching Only, Bleaching Only Phase 2 (No Intervention).
# Column renaming aligns with downstream variable names used in this script.
catch_scenarios_raw <- read_csv(
  "data/inputs/catch/results_table.csv",
  show_col_types = FALSE
)

catch_scenarios <- catch_scenarios_raw %>%
  select(Location, Scenario = scenario, catch_scenario_ton = estimated_catch_tons)

# 2b. Fisher interview data (effort, costs, prices, species allocation)
interviews <- read_csv(
  "data/inputs/surveys/finfish_survey_data_cpue.csv",
  show_col_types = FALSE
)

# ============================================================
# 3. Species composition (p_s)
# Catch_s = Catch × p_s
# Source: semi-structured fisher interviews — relative catch
#         proportions averaged across all interviewed fishers
# ============================================================
species_proportions <- tibble(
  species    = c("snapper", "grouper", "jacks", "barracuda", "snook", "grunt"),
  proportion = c(0.75134,   0.11282,   0.065343, 0.063556,   0.003369, 0.003573)
)
# Verify proportions sum to 1
stopifnot(abs(sum(species_proportions$proportion) - 1) < 1e-4)

# ============================================================
# 4. Social allocation parameters (α_s)
# Sales_s      = Catch_s × α_s
# Consumption_s = Catch_s × (1 − α_s − share_s − waste_s)
# Source: fisher interviews — destination shares by species & location
# ============================================================
destination_shares <- tribble(
  ~Location,      ~species,    ~sale, ~consumption, ~share, ~waste,

  # --- Caye Caulker ---
  "Caye Caulker", "snapper",   0.90,  0.05,         0.05,   0.00,
  "Caye Caulker", "grouper",   1.00,  0.00,         0.00,   0.00,
  "Caye Caulker", "barracuda", 0.90,  0.00,         0.10,   0.00,
  "Caye Caulker", "jacks",     1.00,  0.00,         0.00,   0.00,
  "Caye Caulker", "snook",     0.00,  0.00,         1.00,   0.00,
  "Caye Caulker", "grunt",     0.00,  0.00,         1.00,   0.00,

  # --- Placencia ---
  "Placencia",    "snapper",   0.80,  0.10,         0.05,   0.00,
  "Placencia",    "grouper",   0.80,  0.10,         0.10,   0.00,
  "Placencia",    "barracuda", 0.60,  0.30,         0.00,   0.10,
  "Placencia",    "jacks",     0.35,  0.30,         0.35,   0.00,
  "Placencia",    "snook",     0.35,  0.58,         0.07,   0.00,
  "Placencia",    "grunt",     0.46,  0.30,         0.18,   0.06
)

# ============================================================
# 4b. Species attribution summary & visualization
#     Uses destination_shares as the DEFINITIVE source (defined above).
#     DO NOT redefine destination_shares here.
# ============================================================

# Summary table: contribution of each species × location to each dimension
species_attribution <- destination_shares %>%
  pivot_longer(
    cols      = c(sale, consumption, share, waste),
    names_to  = "dimension",
    values_to = "proportion"
  ) %>%
  filter(proportion > 0) %>%
  mutate(
    dimension = recode(dimension,
      sale        = "Economic (Sale)",
      consumption = "Nutrition (Self-consumption)",
      share       = "Social Cohesion (Community sharing)",
      waste       = "Post-harvest loss"
    ),
    species   = str_to_title(species),
    Location  = factor(Location, levels = c("Caye Caulker", "Placencia"))
  )

# Print attribution table
cat("\n=== SPECIES ATTRIBUTION BY DIMENSION AND LOCALITY ===\n")
print(
  species_attribution %>%
    pivot_wider(names_from = dimension, values_from = proportion, values_fill = 0) %>%
    arrange(Location, species)
)

# Grouped bar plot: species contributions by locality and dimension
p_attribution <- species_attribution %>%
  filter(dimension != "Post-harvest loss") %>%   # keep 3 main dimensions
  ggplot(aes(
    x    = fct_reorder(species, proportion, .desc = TRUE),
    y    = proportion,
    fill = Location
  )) +
  geom_col(position = "dodge", alpha = 0.85) +
  facet_wrap(
    ~ dimension,
    nrow   = 1,
    scales = "free_y",
    labeller = label_wrap_gen(width = 20)
  ) +
  scale_fill_manual(
    values = c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title    = "Species Contributions by Locality and Value Chain Dimension",
    subtitle = "Source: destination_shares (fisher interviews) — proportions of species catch allocated per destination",
    y        = "Proportion of species catch",
    x        = "Species",
    fill     = "Location"
  )

ggsave("data/outputs/fig_species_attribution.png", p_attribution,
       width = 13, height = 5, dpi = 300)

# ============================================================
# 5. Ex-vessel prices (USD/lb) by species & location
# Source: fisher-reported prices from interviews
# ============================================================
average_prices <- interviews %>%
  filter(
    fishing_type_main %in% c("Commercial_Fishing"),
    Location %in% c("Caye Caulker", "Placencia")
  ) %>%
  group_by(Location) %>%
  summarise(
    snapper   = mean(price_snapper_whole, na.rm = TRUE),
    grouper   = mean(price_grouper,       na.rm = TRUE),
    jacks     = mean(price_jacks,         na.rm = TRUE),
    barracuda = mean(price_barracuda,     na.rm = TRUE),
    snook     = mean(price_snook,         na.rm = TRUE),
    grunt     = mean(price_grunt,         na.rm = TRUE),
    .groups   = "drop"
  )

price_long <- average_prices %>%
  pivot_longer(
    cols      = -Location,
    names_to  = "species",
    values_to = "price_lb"
  )

# ============================================================
# 6. Operating costs & annual effort
# Profit = Total_Revenue − (C_trip × N_trips)
# Source: fisher interviews
# ============================================================

# Annual trips per fisher
annual_trips_df <- interviews %>%
  filter(
    fishing_type_main %in% c("Commercial_Fishing"),
    Location %in% c("Caye Caulker", "Placencia"),
    !is.na(number_trips_fish)
  ) %>%
  group_by(Location) %>%
  summarise(
    avg_trips_per_month = mean(number_trips_fish, na.rm = TRUE),
    avg_months_fishing  = mean(months_fishing, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(annual_trips = avg_trips_per_month * avg_months_fishing)

# Operating costs derived cost columns (ice, crew, fuel, food).
# trip cost as: trip_cost = ice + crew + fuel + food 
# trip_revenue = total_catch_lb * weighted_avg_price
# Annual operating cost = trip_cost * annual_trips
#
# Weighted average price per lb (across species, using snapper as proxy
# since it dominates 75% of catch)
avg_price_per_lb <- price_long %>%
  filter(species == "snapper") %>%
  group_by(Location) %>%
  summarise(avg_price = mean(price_lb, na.rm = TRUE), .groups = "drop")

average_costs <- interviews %>%
  filter(fishing_type_main=="Commercial_Fishing",
         Location %in% c("Caye Caulker","Placencia")) %>%
  group_by(Location) %>%
  summarise(
    cost_ice  = mean(cost_ice,   na.rm=TRUE),
    cost_crew = mean(cost_crew, na.rm=TRUE),
    cost_fuel = mean(cost_fuel,  na.rm=TRUE),
    cost_food = mean(cost_Food,  na.rm=TRUE),
    #cost_gear = mean(cost_gear, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(trip_cost = cost_ice + cost_crew + cost_fuel + cost_food) %>%
  left_join(annual_trips_df %>% select(Location, annual_trips), by="Location") %>%
  mutate(annual_operating_cost = trip_cost * annual_trips)



# ============================================================
# 7. Multispecies catch partitioning
# Catch_s = catch_scenario_ton × TON_TO_LB × p_s
# ============================================================
catch_species <- catch_scenarios %>%
  dplyr::mutate(annual_catch_lb = catch_scenario_ton * TON_TO_LB) %>%
  tidyr::crossing(species_proportions) %>%
  dplyr::mutate(
    species_catch_lb = annual_catch_lb * proportion,
    Scenario = factor(Scenario, levels = SCENARIO_ORDER)
  )

# ============================================================
# 8. Social allocation
# Sales_s      = species_catch_lb × α_s
# Consumption_s = species_catch_lb × (1 − α_s)
# ============================================================
economic_base <- catch_species %>%
  left_join(destination_shares, by = c("Location", "species")) %>%
  mutate(
    sale_lb        = species_catch_lb * sale,
    consumption_lb = species_catch_lb * consumption,
    share_lb       = species_catch_lb * share,
    waste_lb       = species_catch_lb * waste
  )

# ============================================================
# 9. Economic valuation
# Revenue = Σ(Sales_s × Price_s)
# ============================================================
economic_results <- economic_base %>%
  left_join(price_long, by = c("Location", "species")) %>%
  mutate(
    revenue           = sale_lb * price_lb,
    consumption_value = consumption_lb * price_lb   # monetary value of self-consumed fish
  )

# Save intermediate output for use by restaurant script
write_csv(
  economic_results,
  "data/outputs/economic_results.csv"
)

# ============================================================
# 10. Aggregate revenue by scenario
# ============================================================
revenue_summary <- economic_results %>%
  group_by(Location, Scenario) %>%
  summarise(
    total_revenue           = sum(revenue,           na.rm = TRUE),
    total_consumption_value = sum(consumption_value, na.rm = TRUE),
    total_sale_lb           = sum(sale_lb,           na.rm = TRUE),
    total_consumption_lb    = sum(consumption_lb,    na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 11. Food security metric
# Meals = Σ(Consumption_s) / LB_PER_MEAL
# ============================================================
consumption_summary <- revenue_summary %>%
  mutate(
    meals_equivalent = total_consumption_lb / LB_PER_MEAL
  )

# ============================================================
# 12. Profit
# Profit = Revenue − (C_trip × N_trips)
# NOTE: Household consumption value is intentionally EXCLUDED
#       from profit. It is reported separately as an in-kind
#       food-security benefit to avoid double-counting.
#       Revenue only accounts for sold catch.
# ============================================================
final_results <- consumption_summary %>%
  left_join(
    average_costs %>% select(Location, annual_operating_cost),
    by = "Location"
  ) %>%
  mutate(
    profit = total_revenue - annual_operating_cost,   # REVISED: no consumption_value
    Scenario = factor(Scenario, levels = SCENARIO_ORDER)
  ) %>%
  arrange(Location, Scenario)

# ── Diagnostic check: is Caye Caulker Baseline profit negative? ──────────────
# This is an important validity check, not just a flag.
# Caye Caulker fishers report very high operating costs (crew, fuel, ice)
# relative to finfish CPUE. If finfish-only profit is negative, this is
# CONSISTENT with interview evidence: most Caye Caulker fishers are
# lobster/conch specialists who catch finfish opportunistically.
# Finfish revenue alone does not cover trip costs — lobster & conch
# subsidise the operation. This result should be REPORTED, not corrected.
cc_baseline_profit <- final_results %>%
  filter(Location == "Caye Caulker", Scenario == "Baseline") %>%
  pull(profit)

if (length(cc_baseline_profit) > 0 && !is.na(cc_baseline_profit) &&
    cc_baseline_profit < 0) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════╗\n")
  cat("║  ⚠️  VALIDITY CHECK — Caye Caulker Baseline Profit is NEGATIVE  ║\n")
  cat("╠══════════════════════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Finfish-only profit: USD %.0f / year                          ║\n",
              cc_baseline_profit))
  cat("║                                                                  ║\n")
  cat("║  This result is EXPECTED and CONSISTENT with field evidence:     ║\n")
  cat("║  • Caye Caulker fishers are primarily lobster/conch specialists  ║\n")
  cat("║  • Finfish is caught opportunistically during lobster trips       ║\n")
  cat("║  • High crew + fuel costs are shared across target species       ║\n")
  cat("║  • Finfish revenue alone does not cover full trip costs          ║\n")
  cat("║                                                                  ║\n")
  cat("║  ACTION: Report this result. Add discussion note in the Rmd     ║\n")
  cat("║  explaining multi-species livelihood complementarity.            ║\n")
  cat("╚══════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
} else if (length(cc_baseline_profit) > 0 && !is.na(cc_baseline_profit) &&
           cc_baseline_profit >= 0) {
  cat("\n✅ Caye Caulker Baseline finfish profit is positive (USD",
      round(cc_baseline_profit, 0), "/ year).\n")
  cat("   Verify: are operating costs correctly attributed to finfish trips only?\n\n")
}

# ============================================================
# 12b. Two-panel temporal table
#      No Intervention : Baseline, Baseline, Bleaching Only, Bleaching Only Phase 2
#      With Restoration: Baseline, Phase 1,  Bleaching,       Phase 2
#      Mirrors the pattern used in catch_estimation.R
# ============================================================
final_results_2panel_temporal <- bind_rows(

  # ---- No Intervention panel (4 synthetic time points) ----
  final_results %>%
    filter(Scenario == "Baseline") %>%
    mutate(Scenario = as.character(Scenario), Scenario = "Year zero",
           scenario_group = "No Intervention"),
  final_results %>%
    filter(Scenario == "Baseline") %>%
    mutate(Scenario = as.character(Scenario), Scenario = "1 Year",
           scenario_group = "No Intervention"),
  final_results %>%
    filter(Scenario == "Bleaching Only") %>%
    mutate(Scenario = as.character(Scenario), Scenario = "2 Years",
           scenario_group = "No Intervention"),
  final_results %>%
    filter(Scenario == "Bleaching Only Phase 2") %>%
    mutate(Scenario = as.character(Scenario), Scenario = "~5 Years",
           scenario_group = "No Intervention"),

  # ---- With Restoration panel (renamed to temporal labels) ----
  final_results %>%
    filter(Scenario %in% c("Baseline", "Phase 1", "Bleaching", "Phase 2")) %>%
    mutate(
      Scenario = recode(as.character(Scenario),
                        "Baseline"  = "Year zero",
                        "Phase 1"   = "1 Year",
                        "Bleaching" = "2 Years",
                        "Phase 2"   = "~5 Years"),
      scenario_group = "With Restoration"
    )

) %>%
  mutate(
    Scenario       = factor(Scenario, levels = c("Year zero", "1 Year", "2 Years", "~5 Years")),
    scenario_group = factor(scenario_group, levels = c("No Intervention", "With Restoration")),
    Location       = factor(Location)
  )

# ============================================================
# 13. Scenario-to-scenario Δ (delta) calculations
#     Computed within each scenario_group × Location so that
#     lag() never crosses panel boundaries.
# ============================================================
delta_results <- final_results_2panel_temporal %>%
  group_by(scenario_group, Location) %>%
  arrange(Scenario, .by_group = TRUE) %>%
  mutate(
    delta_revenue        = total_revenue        - lag(total_revenue),
    delta_profit         = profit               - lag(profit),
    delta_consumption_lb = total_consumption_lb - lag(total_consumption_lb),
    delta_meals          = meals_equivalent     - lag(meals_equivalent),
    transition = case_when(
      Scenario == "1 Year"   ~ "Year zero → 1 Year",
      Scenario == "2 Years"  ~ "1 Year → 2 Years",
      Scenario == "~5 Years" ~ "2 Years → ~5 Years",
      TRUE                   ~ NA_character_
    )
  ) %>%
  ungroup()

# ============================================================
# 14. Year zero → ~5 Years net benefit — both panels
#     With Restoration : Baseline → Phase 2
#     No Intervention  : Baseline → Bleaching Only Phase 2
# ============================================================
restoration_delta <- final_results_2panel_temporal %>%
  filter(Scenario %in% c("Year zero", "~5 Years")) %>%
  group_by(scenario_group, Location) %>%
  summarise(
    baseline_revenue        = total_revenue[Scenario == "Year zero"],
    final_revenue           = total_revenue[Scenario == "~5 Years"],
    delta_revenue           = final_revenue - baseline_revenue,
    pct_revenue_change      = (delta_revenue / baseline_revenue) * 100,

    baseline_profit         = profit[Scenario == "Year zero"],
    final_profit            = profit[Scenario == "~5 Years"],
    delta_profit            = final_profit - baseline_profit,
    pct_profit_change       = (delta_profit / baseline_profit) * 100,

    baseline_consumption_lb = total_consumption_lb[Scenario == "Year zero"],
    final_consumption_lb    = total_consumption_lb[Scenario == "~5 Years"],
    delta_consumption_lb    = final_consumption_lb - baseline_consumption_lb,
    pct_consumption_change  = (delta_consumption_lb / baseline_consumption_lb) * 100,

    baseline_meals          = meals_equivalent[Scenario == "Year zero"],
    final_meals             = meals_equivalent[Scenario == "~5 Years"],
    delta_meals             = final_meals - baseline_meals,
    .groups = "drop"
  )

# ============================================================
# 15. Save outputs
# ============================================================
write_csv(final_results,                 "data/outputs/fisher_economics_by_scenario.csv")
write_csv(final_results_2panel_temporal, "data/outputs/fisher_economics_2panel_temporal.csv")
write_csv(delta_results,                 "data/outputs/fisher_economics_deltas.csv")
write_csv(restoration_delta,             "data/outputs/restoration_net_benefit.csv")

# ============================================================
# 16. Plots
# ============================================================

# Location colour palette
loc_levels <- levels(final_results_2panel_temporal$Location)
loc_cols   <- c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")

# National totals for dashed overlay
national_econ_2panel <- final_results_2panel_temporal %>%
  group_by(scenario_group, Scenario) %>%
  summarise(
    profit           = sum(profit,           na.rm = TRUE),
    meals_equivalent = sum(meals_equivalent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Location = factor("National"))

# 16a. Two-panel profit trajectory
p_profit_2panel <- ggplot(
  final_results_2panel_temporal,
  aes(x = Scenario, y = profit, group = Location, color = Location)
) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +
  geom_line(
    data = national_econ_2panel,
    aes(x = Scenario, y = profit, group = 1, linetype = "National total"),
    linewidth = 1.2, color = "black"
  ) +
  geom_point(
    data = national_econ_2panel,
    aes(x = Scenario, y = profit),
    size = 3, shape = 18, color = "black"
  ) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_color_manual(values = loc_cols) +
  scale_linetype_manual(
    name   = NULL,
    values = c("National total" = "dashed"),
    guide  = guide_legend(
      override.aes = list(color = "black", linewidth = 1.2, shape = NA)
    )
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 20, hjust = 1),
    axis.title         = element_text(face = "bold"),
    legend.position    = "right",
    strip.text         = element_text(face = "bold", size = 13)
  ) +
  labs(
    x     = NULL,
    y     = "Annual Profit (USD)",
    color = "Location"
  )

ggsave("data/outputs/fig_profit_2panel.png", p_profit_2panel,
       width = 12, height = 5.2, dpi = 300)

# 16b. Two-panel meals-equivalent trajectory
p_meals_2panel <- ggplot(
  final_results_2panel_temporal,
  aes(x = Scenario, y = meals_equivalent, group = Location, color = Location)
) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +
  geom_line(
    data = national_econ_2panel,
    aes(x = Scenario, y = meals_equivalent, group = 1, linetype = "National total"),
    linewidth = 1.2, color = "black"
  ) +
  geom_point(
    data = national_econ_2panel,
    aes(x = Scenario, y = meals_equivalent),
    size = 3, shape = 18, color = "black"
  ) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_color_manual(values = loc_cols) +
  scale_linetype_manual(
    name   = NULL,
    values = c("National total" = "dashed"),
    guide  = guide_legend(
      override.aes = list(color = "black", linewidth = 1.2, shape = NA)
    )
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 20, hjust = 1),
    axis.title         = element_text(face = "bold"),
    legend.position    = "right",
    strip.text         = element_text(face = "bold", size = 13)
  ) +
  labs(
    title    = "Household Meals Supported by Self-Consumed Fish",
    subtitle = "1 meal = 0.5 lb (227 g) per adult portion",
    x        = NULL,
    y        = "Annual Meals Equivalent",
    color    = "Location"
  )

ggsave("data/outputs/fig_meals_2panel.png", p_meals_2panel,
       width = 12, height = 5.2, dpi = 300)

# ============================================================
# 17. Print summary
# ============================================================
cat("\n=== FISHER ECONOMICS SUMMARY (all scenarios) ===\n")
print(final_results %>%
        select(Location, Scenario, total_revenue, profit, meals_equivalent))

cat("\n=== FISHER ECONOMICS — TWO-PANEL TEMPORAL VIEW ===\n")
print(final_results_2panel_temporal %>%
        select(scenario_group, Location, Scenario, total_revenue, profit, meals_equivalent))

cat("\n=== NET BENEFIT (Year zero → ~5 Years, by panel) ===\n")
print(restoration_delta %>%
        select(scenario_group, Location, delta_revenue, pct_revenue_change,
               delta_profit, pct_profit_change,
               delta_meals))
