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

# Scenario display order
SCENARIO_ORDER <- c("Baseline", "Phase 1", "Bleaching", "Phase 2")

# ============================================================
# 2. Load data
# ============================================================

# 2a. Catch & biomass by scenario (PRIMARY CATCH INPUT)
catch_scenarios <- read_csv(
  "data/inputs/catch/location_scenario_results_formatted.csv",
  show_col_types = FALSE
)

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
# 13. Scenario-to-scenario Δ (delta) calculations
# ============================================================
delta_results <- final_results %>%
  group_by(Location) %>%
  arrange(Scenario, .by_group = TRUE) %>%
  mutate(
    delta_revenue        = total_revenue        - lag(total_revenue),
    delta_profit         = profit               - lag(profit),
    delta_consumption_lb = total_consumption_lb - lag(total_consumption_lb),
    delta_meals          = meals_equivalent     - lag(meals_equivalent),
    transition = case_when(
      Scenario == "Phase 1"   ~ "Baseline → Phase 1",
      Scenario == "Bleaching" ~ "Phase 1 → Bleaching",
      Scenario == "Phase 2"   ~ "Bleaching → Phase 2",
      TRUE                    ~ NA_character_
    )
  ) %>%
  ungroup()

# ============================================================
# 14. Baseline → Phase 2 net restoration benefit
# ============================================================
restoration_delta <- final_results %>%
  filter(Scenario %in% c("Baseline", "Phase 2")) %>%
  group_by(Location) %>%
  summarise(
    baseline_revenue        = total_revenue[Scenario == "Baseline"],
    phase2_revenue          = total_revenue[Scenario == "Phase 2"],
    delta_revenue           = phase2_revenue - baseline_revenue,
    pct_revenue_change      = (delta_revenue / baseline_revenue) * 100,

    baseline_profit         = profit[Scenario == "Baseline"],
    phase2_profit           = profit[Scenario == "Phase 2"],
    delta_profit            = phase2_profit - baseline_profit,
    pct_profit_change       = (delta_profit / baseline_profit) * 100,

    baseline_consumption_lb = total_consumption_lb[Scenario == "Baseline"],
    phase2_consumption_lb   = total_consumption_lb[Scenario == "Phase 2"],
    delta_consumption_lb    = phase2_consumption_lb - baseline_consumption_lb,
    pct_consumption_change  = (delta_consumption_lb / baseline_consumption_lb) * 100,

    baseline_meals          = meals_equivalent[Scenario == "Baseline"],
    phase2_meals            = meals_equivalent[Scenario == "Phase 2"],
    delta_meals             = phase2_meals - baseline_meals,
    .groups = "drop"
  )

# ============================================================
# 15. Save outputs
# ============================================================
write_csv(final_results,     "data/outputs/fisher_economics_by_scenario.csv")
write_csv(delta_results,     "data/outputs/fisher_economics_deltas.csv")
write_csv(restoration_delta, "data/outputs/restoration_net_benefit.csv")

# ============================================================
# 16. Plots
# ============================================================

# 16a. Profit trajectory across scenarios
p_profit <- ggplot(
  final_results,
  aes(x = Scenario, y = profit, group = Location, color = Location)
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")) +
  theme_classic(base_size = 13) +
  labs(
    title    = "Annual Fisher Profit Under Coral Restoration Scenarios",
    subtitle = "Profit = Revenue − Variable Operating Costs",
    y        = "Annual Profit (BZD)",
    x        = NULL,
    color    = "Location"
  )

ggsave("data/outputs/fig_profit_trajectory.png", p_profit,
       width = 8, height = 5, dpi = 300)

# 16b. Meals equivalent across scenarios
p_meals <- ggplot(
  final_results,
  aes(x = Scenario, y = meals_equivalent, group = Location, color = Location)
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")) +
  theme_classic(base_size = 13) +
  labs(
    title    = "Household Meals Supported by Self-Consumed Fish",
    subtitle = "1 meal = 0.5 lb (227 g) per adult portion",
    y        = "Annual Meals Equivalent",
    x        = NULL,
    color    = "Location"
  )

ggsave("data/outputs/fig_meals_trajectory.png", p_meals,
       width = 8, height = 5, dpi = 300)

# ============================================================
# 17. Print summary
# ============================================================
cat("\n=== FISHER ECONOMICS SUMMARY ===\n")
print(final_results %>%
        select(Location, Scenario, total_revenue, profit, meals_equivalent))

cat("\n=== NET RESTORATION BENEFIT (Baseline → Phase 2) ===\n")
print(restoration_delta %>%
        select(Location, delta_revenue, pct_revenue_change,
               delta_profit, pct_profit_change,
               delta_meals))
