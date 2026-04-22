# ============================================================
# Script 02: Restaurant Seafood Demand
# Project:   Coral Restoration Benefits — Belize SSF
# Purpose:   Aggregate annual restaurant seafood demand by
#            species and location; compare with fisher supply;
#            price discrepancy analysis; per-species exploratory
# Author:    [Your name]
# Date:      2026
# ============================================================

# ============================================================
# 0. Packages
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)
library(forcats)
library(scales)

# ============================================================
# 1. Constants
# ============================================================
HIGH_MONTHS    <- 6   # December–May (high season)
LOW_MONTHS     <- 6   # June–November (low season)
DAYS_PER_MONTH <- 30
SCENARIO_ORDER <- c("Baseline", "Phase 1", "Bleaching", "Phase 2",
                    "Bleaching Only", "Bleaching Only Phase 2")

# ============================================================
# 2. Load data
# ============================================================
df_raw <- read_excel(
  "data/inputs/surveys/restaurants_answers_survey123.xlsx",
  sheet = "survey",
  na    = "NA"
)

interviews <- read_csv(
  "data/inputs/surveys/finfish_survey_data_cpue.csv",
  show_col_types = FALSE
)

economic_results <- read_csv(
  "data/outputs/economic_results.csv",
  show_col_types = FALSE
)

# ============================================================
# 3. Sampling design note
# ============================================================
# IMPORTANT: The two localities have fundamentally different
# sampling designs and this MUST be acknowledged before any
# comparison of absolute volumes.
#
# Caye Caulker: n = 16 restaurants surveyed.
#   The total number of restaurants in Caye Caulker is not
#   precisely known. The sample likely represents a partial
#   coverage of the restaurant sector. Therefore:
#   - Absolute demand totals for Caye Caulker are LOWER BOUNDS,
#     not population estimates.
#   - We do NOT extrapolate to a population total because we
#     lack a reliable denominator.
#   - Comparisons with Placencia should focus on per-restaurant
#     averages (mean_demand_per_restaurant), not totals.
#
# Placencia: n = 31 restaurants surveyed.
#   This is understood to represent a near-complete census of
#   the restaurant sector in Placencia village.
#   - Absolute demand totals are treated as population-level.
#
# Analytical approach:
#   - Report both total and per-restaurant demand.
#   - Flag Caye Caulker totals as partial-sample estimates.
#   - Use per-restaurant averages for cross-locality comparisons.

sampling_note <- tibble(
  Location          = c("Caye Caulker", "Placencia"),
  n_surveyed        = c(16L, 31L),
  sample_type       = c("Partial sample — total unknown",
                        "Near-complete census"),
  use_total_demand  = c(FALSE, TRUE),
  comparison_metric = c("Per-restaurant average only",
                        "Total AND per-restaurant average")
)

cat("\n=== SAMPLING DESIGN CONTEXT ===\n")
print(sampling_note)
cat("\n⚠️  Caye Caulker: n=16 restaurants. Total population unknown.\n")
cat("   Do NOT compare absolute demand volumes between localities.\n")
cat("   Use per-restaurant averages for cross-locality comparisons.\n\n")

# ============================================================
# 4. Clean restaurant data
# ============================================================
parse_num <- function(x) readr::parse_number(as.character(x))

df_clean <- df_raw %>%
  filter(Site %in% c("Caye Caulker", "Placencia")) %>%
  mutate(
    Site                      = factor(Site, levels = c("Caye Caulker", "Placencia")),
    customers_high            = as.numeric(customers_high),
    customers_low             = as.numeric(customers_low),
    months_operate_restaurant = as.numeric(months_operate_restaurant),
    across(starts_with("buy_"),   parse_num),
    across(starts_with("price_"), parse_num),
    # Selling prices (what restaurants charge customers)
    across(starts_with("selling_"), parse_num),
    portion_size_Finfish = parse_num(portion_size_Finfish)
  )

# ============================================================
# 5. Annual customer volume
# ============================================================
df_clean <- df_clean %>%
  mutate(
    annual_customers =
      (customers_high * HIGH_MONTHS * DAYS_PER_MONTH +
       customers_low  * LOW_MONTHS  * DAYS_PER_MONTH) *
      (months_operate_restaurant / 12)
  )

# ============================================================
# 6. Purchase frequency multiplier
# ============================================================
df_clean <- df_clean %>%
  mutate(
    purchase_multiplier = case_when(
      frequency_purchase == "Daily"   ~ 365,
      frequency_purchase == "Weekly"  ~ 52,
      frequency_purchase == "Monthly" ~ 12,
      TRUE                            ~ NA_real_
    )
  )

# ============================================================
# 7. Annual seafood demand by species
# ============================================================
df_clean <- df_clean %>%
  mutate(
    annual_snapper   =
      coalesce(buy_Snapper_low,          0) * purchase_multiplier * (LOW_MONTHS  / 12) +
      coalesce(buy_Snapper_high_season,  0) * purchase_multiplier * (HIGH_MONTHS / 12),
    annual_grouper   =
      coalesce(buy_Grouper_low,          0) * purchase_multiplier * (LOW_MONTHS  / 12) +
      coalesce(buy_Grouper_high_season,  0) * purchase_multiplier * (HIGH_MONTHS / 12),
    annual_snook     =
      coalesce(buy_Snook_low,            0) * purchase_multiplier * (LOW_MONTHS  / 12) +
      coalesce(buy_Snook_high_season,    0) * purchase_multiplier * (HIGH_MONTHS / 12),
    annual_barracuda =
      coalesce(buy_Barracuda_low,        0) * purchase_multiplier * (LOW_MONTHS  / 12) +
      coalesce(buy_Barracuda_high_season,0) * purchase_multiplier * (HIGH_MONTHS / 12),
    annual_jacks     =
      coalesce(buy_Jacks_low,            0) * purchase_multiplier * (LOW_MONTHS  / 12) +
      coalesce(buy_Jacks_high_season,    0) * purchase_multiplier * (HIGH_MONTHS / 12),
    total_annual_seafood =
      annual_snapper + annual_grouper + annual_snook +
      annual_barracuda + annual_jacks
  )

# ============================================================
# 8. Per-species exploratory analysis (BEFORE aggregation)
#    Quantity, portion size, purchase price, selling price,
#    and implied markup per species per restaurant
# ============================================================

# 8a. Per-restaurant, per-species demand breakdown
species_demand_long <- df_clean %>%
  select(Site, name_estaurant, type_restaurant,
         annual_snapper, annual_grouper, annual_snook,
         annual_barracuda, annual_jacks,
         price_Snapper, price_Grouper, price_Snook,
         price_Barracuda, price_Jacks,
         selling_Snapper, selling_Grouper, selling_Jacks,
         selling_Barracuda,
         portion_size_Finfish) %>%
  pivot_longer(
    cols      = starts_with("annual_"),
    names_to  = "species",
    values_to = "annual_qty_lb"
  ) %>%
  mutate(
    species = str_remove(species, "annual_") %>% str_to_title(),
    # Match purchase price column to species
    purchase_price_lb = case_when(
      species == "Snapper"   ~ price_Snapper,
      species == "Grouper"   ~ price_Grouper,
      species == "Snook"     ~ price_Snook,
      species == "Barracuda" ~ price_Barracuda,
      species == "Jacks"     ~ price_Jacks,
      TRUE                   ~ NA_real_
    ),
    # Match selling price column to species
    # NOTE: selling_Snook is absent from the survey instrument —
    #       restaurants were not asked the selling price for snook dishes.
    #       We assign NA for snook selling price.
    selling_price_dish = case_when(
      species == "Snapper"   ~ selling_Snapper,
      species == "Grouper"   ~ selling_Grouper,
      species == "Snook"     ~ NA_real_,
      species == "Barracuda" ~ selling_Barracuda,
      species == "Jacks"     ~ selling_Jacks,
      TRUE                   ~ NA_real_
    ),
    # Portion size in lb (convert from oz if needed — stored as numeric oz)
    portion_lb = portion_size_Finfish / 16,
    # Cost of fish per dish = purchase_price_lb × portion_lb
    fish_cost_per_dish = purchase_price_lb * portion_lb,
    # Markup = (selling_price - fish_cost) / fish_cost
    # Note: this is a PARTIAL markup — other costs (labour, overhead) not included
    markup_ratio = (selling_price_dish - fish_cost_per_dish) / fish_cost_per_dish
  ) %>%
  filter(!is.na(annual_qty_lb), annual_qty_lb > 0)

# 8b. Summary by species and location
species_demand_summary <- species_demand_long %>%
  group_by(Site, species) %>%
  summarise(
    n_restaurants        = n(),
    mean_annual_qty_lb   = mean(annual_qty_lb,      na.rm = TRUE),
    total_annual_qty_lb  = sum(annual_qty_lb,        na.rm = TRUE),
    mean_purchase_price  = mean(purchase_price_lb,   na.rm = TRUE),
    mean_selling_price   = mean(selling_price_dish,  na.rm = TRUE),
    mean_portion_oz      = mean(portion_size_Finfish,na.rm = TRUE),
    mean_markup_ratio    = mean(markup_ratio,         na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== PER-SPECIES RESTAURANT DEMAND EXPLORATORY SUMMARY ===\n")
print(species_demand_summary)

write_csv(species_demand_summary, "data/outputs/restaurant_species_demand_exploratory.csv")
write_csv(species_demand_long,    "data/outputs/restaurant_species_demand_per_restaurant.csv")

# Plot: annual demand by species and location (per-restaurant average)
p_species_per_rest <- ggplot(
  species_demand_summary,
  aes(x = fct_reorder(species, mean_annual_qty_lb, .desc = TRUE),
      y = mean_annual_qty_lb, fill = Site)
) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")) +
  scale_y_continuous(labels = comma) +
  theme_classic(base_size = 12) +
  labs(
    title    = "Mean Annual Seafood Demand per Restaurant by Species",
    subtitle = "Per-restaurant average — comparable across localities regardless of sample size",
    y        = "Mean annual demand (lb/restaurant/year)",
    x        = "Species",
    fill     = "Location"
  )

ggsave("data/outputs/fig_species_demand_per_restaurant.png", p_species_per_rest,
       width = 9, height = 5, dpi = 300)

# ============================================================
# 9. Price discrepancy analysis
#    Fisher ex-vessel price vs. restaurant purchase price
#
#    IMPORTANT CONTEXT: The column source_fish in the restaurant
#    survey shows that the vast majority of restaurants purchase
#    directly from local fishers (no intermediary). Therefore,
#    any price gap is NOT a middleman margin — it reflects:
#      (a) Negotiation dynamics in direct fisher-restaurant sales
#      (b) Quality/form differences (whole fish vs. fillet)
#      (c) Seasonal price variation
#      (d) Possible recall bias in either survey
#    We explore the gap descriptively to understand these dynamics.
# ============================================================

# Fisher ex-vessel prices (from interviews)
fisher_prices <- interviews %>%
  filter(
    fishing_type_main %in% c("Commercial_Fishing", "Sport_Deep_Sea_Fishing",
                              "Sport_Fly_Fishing"),
    Location %in% c("Caye Caulker", "Placencia")
  ) %>%
  group_by(Location) %>%
  summarise(
    fisher_snapper   = mean(price_snapper_whole, na.rm = TRUE),
    fisher_grouper   = mean(price_grouper,       na.rm = TRUE),
    fisher_snook     = mean(price_snook,         na.rm = TRUE),
    fisher_barracuda = mean(price_barracuda,     na.rm = TRUE),
    fisher_jacks     = mean(price_jacks,         na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-Location, names_to = "species", values_to = "fisher_price_lb") %>%
  mutate(species = str_remove(species, "fisher_") %>% str_to_title())

# Restaurant purchase prices (what restaurants pay)
restaurant_prices <- df_clean %>%
  group_by(Site) %>%
  summarise(
    Snapper   = mean(price_Snapper,   na.rm = TRUE),
    Grouper   = mean(price_Grouper,   na.rm = TRUE),
    Snook     = mean(price_Snook,     na.rm = TRUE),
    Barracuda = mean(price_Barracuda, na.rm = TRUE),
    Jacks     = mean(price_Jacks,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-Site, names_to = "species", values_to = "restaurant_price_lb")

# Merge and compute gap
price_comparison <- fisher_prices %>%
  left_join(restaurant_prices, by = c("Location" = "Site", "species")) %>%
  mutate(
    price_gap_abs = restaurant_price_lb - fisher_price_lb,
    price_gap_pct = (price_gap_abs / fisher_price_lb) * 100,
    gap_direction = case_when(
      price_gap_abs > 0  ~ "Restaurant pays MORE than fisher reports",
      price_gap_abs < 0  ~ "Restaurant pays LESS than fisher reports",
      price_gap_abs == 0 ~ "Prices match",
      TRUE               ~ "Insufficient data"
    )
  ) %>%
  filter(!is.na(fisher_price_lb), !is.na(restaurant_price_lb))

cat("\n=== PRICE DISCREPANCY: FISHER vs. RESTAURANT PURCHASE PRICE ===\n")
cat("Context: Most restaurants buy directly from fishers (source_fish = 'Local fisher').\n")
cat("Any gap reflects negotiation, form differences, or recall bias — not intermediary margins.\n\n")
print(price_comparison %>%
        select(Location, species, fisher_price_lb, restaurant_price_lb,
               price_gap_abs, price_gap_pct, gap_direction))

write_csv(price_comparison, "data/outputs/price_comparison_fisher_restaurant.csv")

# Lollipop chart: price gap by species and location
p_price_gap <- price_comparison %>%
  ggplot(aes(
    x     = fct_reorder(species, price_gap_pct),
    y     = price_gap_pct,
    color = price_gap_pct > 0
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_segment(aes(xend = species, y = 0, yend = price_gap_pct),
               linewidth = 1.2) +
  geom_point(size = 4) +
  coord_flip() +
  facet_wrap(~Location) +
  scale_color_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "#D73027"),
    labels = c("TRUE" = "Restaurant pays more", "FALSE" = "Restaurant pays less"),
    name   = NULL
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top") +
  labs(
    title    = "Price Discrepancy: Fisher Ex-Vessel vs. Restaurant Purchase Price",
    subtitle = "Direct sales (no intermediary) — gap reflects negotiation, form, or recall differences",
    y        = "Price gap (% relative to fisher-reported price)",
    x        = "Species"
  )

ggsave("data/outputs/fig_price_discrepancy.png", p_price_gap,
       width = 10, height = 5, dpi = 300)

# ============================================================
# 10. Aggregate restaurant demand by location
#     (existing logic — unchanged)
# ============================================================
restaurant_demand <- df_clean %>%
  group_by(Site) %>%
  summarise(
    total_demand_lb            = sum(total_annual_seafood, na.rm = TRUE),
    mean_demand_per_restaurant = mean(total_annual_seafood, na.rm = TRUE),
    n_restaurants              = n(),
    total_employees            = sum(employees, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(sampling_note %>% select(Location, sample_type, use_total_demand),
            by = c("Site" = "Location")) %>%
  mutate(
    demand_note = if_else(
      use_total_demand,
      "Census — total is population estimate",
      "Partial sample — total is lower bound only"
    )
  )

# Species-level demand by location (totals)
restaurant_demand_species <- df_clean %>%
  group_by(Site) %>%
  summarise(
    demand_snapper   = sum(annual_snapper,   na.rm = TRUE),
    demand_grouper   = sum(annual_grouper,   na.rm = TRUE),
    demand_snook     = sum(annual_snook,     na.rm = TRUE),
    demand_barracuda = sum(annual_barracuda, na.rm = TRUE),
    demand_jacks     = sum(annual_jacks,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = starts_with("demand_"),
    names_to  = "species",
    values_to = "demand_lb"
  ) %>%
  mutate(species = str_remove(species, "demand_"))

# ============================================================
# 11. Fisher supply by scenario
# ============================================================
fisher_supply <- economic_results %>%
  group_by(Location, Scenario) %>%
  summarise(
    total_sale_lb = sum(sale_lb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Scenario = factor(Scenario, levels = SCENARIO_ORDER))

# ============================================================
# 12. Supply vs. demand comparison
# ============================================================
supply_demand <- fisher_supply %>%
  left_join(restaurant_demand, by = c("Location" = "Site")) %>%
  mutate(
    surplus_deficit_lb = total_sale_lb - total_demand_lb,
    supply_ratio       = total_sale_lb / total_demand_lb,
    n_fishers          = case_when(
      Location == "Caye Caulker" ~ 13,
      Location == "Placencia"    ~ 9
    ),
    sales_per_fisher   = total_sale_lb / n_fishers,
    fishers_required   = total_demand_lb / sales_per_fisher,
    surplus_interpretation = case_when(
      surplus_deficit_lb > 0 ~ "Surplus — extra fish available for community/new markets",
      surplus_deficit_lb < 0 ~ "Deficit — demand exceeds local fisher supply",
      TRUE                   ~ "Balanced"
    )
  )

# ============================================================
# 13. Restoration effect on fishers required
# ============================================================
restoration_fishers <- supply_demand %>%
  filter(Scenario %in% c("Baseline", "Phase 2")) %>%
  select(Location, Scenario, fishers_required, surplus_deficit_lb) %>%
  pivot_wider(
    names_from  = Scenario,
    values_from = c(fishers_required, surplus_deficit_lb)
  ) %>%
  rename(
    fishers_required_baseline = `fishers_required_Baseline`,
    fishers_required_phase2   = `fishers_required_Phase 2`,
    surplus_baseline          = `surplus_deficit_lb_Baseline`,
    surplus_phase2            = `surplus_deficit_lb_Phase 2`
  ) %>%
  mutate(
    reduction_fishers_needed = fishers_required_baseline - fishers_required_phase2,
    delta_surplus_lb         = surplus_phase2 - surplus_baseline
  )

# ============================================================
# 13b. Two-panel temporal supply-demand table
#      No Intervention : Baseline, Baseline, Bleaching Only, Bleaching Only Phase 2
#      With Restoration: Baseline, Phase 1,  Bleaching,      Phase 2
#      Mirrors the pattern used in 01_fisher_economics.R §12b
# ============================================================
supply_demand_2panel <- bind_rows(

  # ---- No Intervention panel ----
  supply_demand %>%
    filter(Scenario == "Baseline") %>%
    mutate(Scenario = "Year zero", scenario_group = "No Intervention"),
  supply_demand %>%
    filter(Scenario == "Baseline") %>%
    mutate(Scenario = "1 Year", scenario_group = "No Intervention"),
  supply_demand %>%
    filter(Scenario == "Bleaching Only") %>%
    mutate(Scenario = "2 Years", scenario_group = "No Intervention"),
  supply_demand %>%
    filter(Scenario == "Bleaching Only Phase 2") %>%
    mutate(Scenario = "~5 Years", scenario_group = "No Intervention"),

  # ---- With Restoration panel ----
  supply_demand %>%
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
# 13c. Delta calculations within each panel
#      Change in surplus/deficit, supply_ratio, fishers_required
#      at each time step (lag() never crosses panel boundaries)
# ============================================================
supply_demand_deltas <- supply_demand_2panel %>%
  group_by(scenario_group, Location) %>%
  arrange(Scenario, .by_group = TRUE) %>%
  mutate(
    delta_surplus_lb     = surplus_deficit_lb - lag(surplus_deficit_lb),
    delta_supply_ratio   = supply_ratio       - lag(supply_ratio),
    delta_fishers_req    = fishers_required   - lag(fishers_required),
    transition = case_when(
      Scenario == "1 Year"   ~ "Year zero → 1 Year",
      Scenario == "2 Years"  ~ "1 Year → 2 Years",
      Scenario == "~5 Years" ~ "2 Years → ~5 Years",
      TRUE                   ~ NA_character_
    )
  ) %>%
  ungroup()

# ============================================================
# 13d. Year zero → ~5 Years net benefit comparison (both panels)
#      With Restoration : Baseline → Phase 2
#      No Intervention  : Baseline → Bleaching Only Phase 2
# ============================================================
supply_demand_net_benefit <- supply_demand_2panel %>%
  filter(Scenario %in% c("Year zero", "~5 Years")) %>%
  group_by(scenario_group, Location) %>%
  summarise(
    baseline_supply_lb     = total_sale_lb[Scenario == "Year zero"],
    final_supply_lb        = total_sale_lb[Scenario == "~5 Years"],
    delta_supply_lb        = final_supply_lb - baseline_supply_lb,
    pct_supply_change      = (delta_supply_lb / baseline_supply_lb) * 100,

    baseline_surplus_lb    = surplus_deficit_lb[Scenario == "Year zero"],
    final_surplus_lb       = surplus_deficit_lb[Scenario == "~5 Years"],
    delta_surplus_lb       = final_surplus_lb - baseline_surplus_lb,

    baseline_supply_ratio  = supply_ratio[Scenario == "Year zero"],
    final_supply_ratio     = supply_ratio[Scenario == "~5 Years"],
    delta_supply_ratio     = final_supply_ratio - baseline_supply_ratio,

    baseline_fishers_req   = fishers_required[Scenario == "Year zero"],
    final_fishers_req      = fishers_required[Scenario == "~5 Years"],
    delta_fishers_req      = final_fishers_req - baseline_fishers_req,
    .groups = "drop"
  )

# ============================================================
# 13e. Freed-up catch from reduced fisher requirement
#      Hypothesis: restoration requires fewer fishers to cover
#      restaurant demand. The catch of those freed-up fishers
#      is no longer committed to restaurants and is instead
#      available for self-consumption, community sharing,
#      or sales to other markets.
#
#      freed_catch = freed_fishers × avg_total_catch_per_fisher
#
#      avg_total_catch_per_fisher uses the *final-scenario* catch
#      (Phase 2 for With Restoration; Bleaching Only Phase 2 for
#      No Intervention), reflecting actual productivity at ~5 Years.
#      Total catch includes all destination channels:
#      sale + self-consumption + sharing + waste.
# ============================================================

# Step 1 — total catch per Location/Scenario (sum across all species
#          and all destination channels = species_catch_lb)
total_catch_by_scenario <- economic_results %>%
  group_by(Location, Scenario) %>%
  summarise(
    total_catch_lb = sum(sale_lb + consumption_lb + share_lb + waste_lb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_fishers = case_when(
      Location == "Caye Caulker" ~ 13,
      Location == "Placencia"    ~ 9
    ),
    avg_catch_per_fisher_lb = total_catch_lb / n_fishers
  )

# Step 2 — map each panel's "~5 Years" endpoint to its source scenario
#           so the correct catch-per-fisher is used
final_scenario_lookup <- tibble(
  scenario_group = c("With Restoration", "No Intervention"),
  Scenario_src   = c("Phase 2",          "Bleaching Only Phase 2")
)

# Step 3 — join and compute freed catch
freed_catch_analysis <- supply_demand_net_benefit %>%
  left_join(final_scenario_lookup, by = "scenario_group") %>%
  left_join(
    total_catch_by_scenario %>%
      select(Location, Scenario, avg_catch_per_fisher_lb),
    by = c("Location", "Scenario_src" = "Scenario")
  ) %>%
  mutate(
    # delta_fishers_req < 0 means fewer fishers needed (restoration benefit)
    # freed_fishers is positive when restoration reduces the required count
    freed_fishers    = -delta_fishers_req,
    freed_catch_lb   = freed_fishers * avg_catch_per_fisher_lb,
    freed_catch_ton  = freed_catch_lb / 2204.62,
    freed_catch_note = case_when(
      freed_fishers > 0 ~
        "Catch available for self-consumption, community sharing, or other markets",
      freed_fishers < 0 ~
        "More fishers required — additional catch committed to restaurant supply",
      TRUE ~ "No change in fisher requirement"
    )
  ) %>%
  select(-Scenario_src)

# ============================================================
# 14. Save outputs
# ============================================================
write_csv(restaurant_demand,          "data/outputs/restaurant_demand_summary.csv")
write_csv(restaurant_demand_species,  "data/outputs/restaurant_demand_by_species.csv")
write_csv(supply_demand,              "data/outputs/supply_demand_comparison.csv")
write_csv(restoration_fishers,        "data/outputs/restoration_fishers_required.csv")
write_csv(supply_demand_2panel,       "data/outputs/supply_demand_2panel_temporal.csv")
write_csv(supply_demand_deltas,       "data/outputs/supply_demand_deltas.csv")
write_csv(supply_demand_net_benefit,  "data/outputs/supply_demand_net_benefit.csv")
write_csv(freed_catch_analysis,       "data/outputs/freed_catch_analysis.csv")

# ============================================================
# 15. Plots
# ============================================================

# Colour palette (consistent with other scripts)
loc_cols <- c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")

# 15a. Restaurant demand by species (total — with sampling caveat label)
p_species_demand <- ggplot(
  restaurant_demand_species,
  aes(x = fct_reorder(species, demand_lb, .desc = TRUE),
      y = demand_lb, fill = Site)
) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Caye Caulker" = "#1B7837", "Placencia" = "#762A83")) +
  scale_y_continuous(labels = comma) +
  theme_classic(base_size = 13) +
  labs(
    title    = "Annual Restaurant Seafood Demand by Species (Total)",
    subtitle = "⚠️ Caye Caulker: partial sample (n=16) — total is lower bound. Placencia: census (n=31)",
    y        = "Annual Demand (lb)",
    x        = "Species",
    fill     = "Location"
  )

ggsave("data/outputs/fig_restaurant_demand_species.png", p_species_demand,
       width = 8, height = 5, dpi = 300)

# 15c. Two-panel supply ratio trajectory
#      (supply_ratio = fisher_supply_lb / restaurant_demand_lb)
#      Values > 1 = fishers can fully cover demand; < 1 = deficit
p_supply_ratio_2panel <- ggplot(
  supply_demand_2panel,
  aes(x = Scenario, y = supply_ratio, group = Location, color = Location)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_color_manual(values = loc_cols) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
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
    title    = "Fisher Supply Capacity Relative to Restaurant Demand",
    subtitle = "Ratio > 1 = supply covers demand; < 1 = deficit\n⚠️ Caye Caulker demand is a partial-sample lower bound",
    x        = NULL,
    y        = "Supply / Demand ratio",
    color    = "Location"
  )

ggsave("data/outputs/fig_supply_ratio_2panel.png", p_supply_ratio_2panel,
       width = 12, height = 5.2, dpi = 300)

# 15d. Two-panel surplus/deficit trajectory (lb/year)
p_surplus_2panel <- ggplot(
  supply_demand_2panel,
  aes(x = Scenario, y = surplus_deficit_lb, group = Location, color = Location)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.6) +
  facet_wrap(~ scenario_group, ncol = 2, scales = "free_x") +
  scale_color_manual(values = loc_cols) +
  scale_y_continuous(labels = comma) +
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
    title    = "Fisher Supply Surplus / Deficit vs. Restaurant Demand",
    subtitle = "Positive = surplus; Negative = demand exceeds fisher supply\n⚠️ Caye Caulker demand is a partial-sample lower bound",
    x        = NULL,
    y        = "Surplus / Deficit (lb/year)",
    color    = "Location"
  )

ggsave("data/outputs/fig_surplus_2panel.png", p_surplus_2panel,
       width = 12, height = 5.2, dpi = 300)

# 15e. Freed-up catch comparison between panels (Year zero → ~5 Years)
#      Bars show catch (lb/year) of freed-up fishers per location and panel
p_freed_catch <- freed_catch_analysis %>%
  filter(freed_fishers > 0) %>%
  ggplot(aes(x = Location, y = freed_catch_lb, fill = scenario_group)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(
    values = c("No Intervention" = "#D73027", "With Restoration" = "#1B7837"),
    name   = NULL
  ) +
  scale_y_continuous(labels = comma) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  labs(
    title    = "Catch Available for Alternative Uses (Year zero → ~5 Years)",
    subtitle = "Catch of fishers freed from restaurant supply obligation\nAvailable for self-consumption, community sharing, or other markets",
    y        = "Freed catch (lb/year)",
    x        = NULL
  )

ggsave("data/outputs/fig_freed_catch.png", p_freed_catch,
       width = 7, height = 5, dpi = 300)

# ============================================================
# 16. Print summary
# ============================================================
cat("\n=== RESTAURANT DEMAND SUMMARY ===\n")
print(restaurant_demand %>%
        select(Site, n_restaurants, total_demand_lb,
               mean_demand_per_restaurant, demand_note))

cat("\n=== SUPPLY vs DEMAND BY SCENARIO ===\n")
print(supply_demand %>%
        select(Location, Scenario, total_sale_lb, total_demand_lb,
               surplus_deficit_lb, supply_ratio))

cat("\n=== RESTORATION EFFECT ON FISHERS REQUIRED ===\n")
print(restoration_fishers)

cat("\n=== SUPPLY vs DEMAND — TWO-PANEL TEMPORAL VIEW ===\n")
print(supply_demand_2panel %>%
        select(scenario_group, Location, Scenario,
               total_sale_lb, surplus_deficit_lb, supply_ratio, fishers_required))

cat("\n=== NET BENEFIT (Year zero → ~5 Years, by panel) ===\n")
print(supply_demand_net_benefit %>%
        select(scenario_group, Location,
               delta_supply_lb, pct_supply_change,
               delta_surplus_lb, delta_supply_ratio, delta_fishers_req))

cat("\n=== FREED-UP CATCH ANALYSIS ===\n")
cat("Fewer fishers needed to cover restaurant demand → their catch available for other uses\n")
print(freed_catch_analysis %>%
        select(scenario_group, Location,
               baseline_fishers_req, final_fishers_req, freed_fishers,
               avg_catch_per_fisher_lb, freed_catch_lb, freed_catch_ton,
               freed_catch_note))
