# Robust Biomass Energy Analysis with Country-Specific Energy Shares
#   • Correctly processes country-level fire impacts
#   • Applies country-specific energy shares
#   • Features enhanced visualizations
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(terra)
library(sf)
library(rnaturalearth)
library(exactextractr)
library(readxl)
library(scales)
library(patchwork)
library(ggtext)
library(glue)
library(ggridges)
library(viridis)
library(ggrepel)

# Configuration ───────────────────────────────────────────────────────────────
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
output_dir <- "H:/UPTAKE/Biomass_Loss_Results"
g4m_dir <- "H:/Uptake"
harvest_file <- "H:/UPTAKE/EU_harvest_residues_selected_scenarios.csv"
energy_share_file <- "H:/Uptake/Energy_share_country.xlsx"
results_dir <- file.path(output_dir, "Fire_Impact_Analysis")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# EU27+UK countries
eu27_uk_countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark",
                       "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
                       "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
                       "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "United Kingdom")

# Biomass conversion factors
wood_density <- 0.5     # t/m³ (oven-dry wood)
carbon_fraction <- 0.5  # fraction of carbon in dry wood
conversion_factor <- wood_density * carbon_fraction  # tC/m³

# Expected value ranges for sanity checks
expected_ranges <- list(
  harvest_m3 = c(min = 1e6, max = 5e8),
  biomass_tc = c(min = 1e6, max = 1e10),
  energy_tc = c(min = 1e5, max = 1e9),
  reduction_pct = c(min = 0, max = 100)
)

# 1. Create EU27+UK boundary and country-level boundaries ────────────────────
countries <- ne_countries(scale = "medium", returnclass = "sf")
eu27_uk_sf <- countries %>% 
  filter(sovereignt %in% eu27_uk_countries) %>%
  st_union() %>% 
  st_sf()

# Country-level boundaries
country_sf <- countries %>% 
  filter(sovereignt %in% eu27_uk_countries) %>%
  select(sovereignt) %>%
  rename(country = sovereignt) %>%
  st_sf()

# 2. Load and process energy share data ───────────────────────────────────────
energy_share_data <- read_excel(energy_share_file) %>%
  rename(Energy_Share = Share_bioenergy)

# Create country code to name mapping
country_mapping <- data.frame(
  Country_Code = c("AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", 
                   "DE", "GR", "HU", "IE", "IT", "LV", "LT", "LU", "MT", "NL", 
                   "PL", "PT", "RO", "SK", "SI", "ES", "SE", "GB"),
  Country = c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark",
              "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
              "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
              "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "United Kingdom")
)

# Create country-specific energy share mapping
country_energy_shares <- energy_share_data %>%
  rename(Country_Code = Country) %>%
  left_join(country_mapping, by = "Country_Code") %>%
  filter(Country %in% eu27_uk_countries) %>%
  select(Country, Energy_Share) %>%
  rename(country = Country, energy_share = Energy_Share)

# 3. CORRECTED: Calculate fire impact with country-specific shares ─────────────
calculate_fire_impact <- function(scenario) {
  # Load harvest data for this scenario
  harvest_scenario <- read_csv(harvest_file) %>%
    filter(Scenario == paste0("Ref_", scenario, "__RCPref__NONE")) %>%
    select(Year, Total_Harvest) %>%
    rename(year = Year, harvest_m3 = Total_Harvest)
  
  # Initialize results with proper column names
  results <- tibble(
    scenario = character(),
    year = integer(),
    total_biomass_tc = double(),
    fire_loss_tc = double(),
    harvest_before_m3 = double(),
    harvest_after_m3 = double(),
    energy_before_tc = double(),
    energy_after_tc = double(),
    reduction_pct = double(),
    prop_loss = double()
  )
  
  for (year in 2015:2100) {
    # Load fire loss raster
    loss_file <- file.path(output_dir, scenario, glue("Total_Loss_{year}.tif"))
    if (!file.exists(loss_file)) next
    
    # Load biomass density raster
    biom_file <- file.path(g4m_dir, scenario, glue("{scenario}_Biom_fm_total_{year}.tif"))
    if (!file.exists(biom_file)) next
    
    tryCatch({
      # Load rasters
      loss_raster <- rast(loss_file)
      biom_raster <- rast(biom_file)
      
      # Create template raster from loss raster
      template <- rast(ext(loss_raster), resolution = res(loss_raster), crs = crs(loss_raster))
      
      # Resample biomass raster to match loss raster
      biom_resampled <- resample(biom_raster, template, method = "near")
      
      # Ensure both rasters have same extent
      loss_cropped <- crop(loss_raster, template)
      biom_cropped <- crop(biom_resampled, template)
      
      # Extract country-level loss and biomass
      country_loss <- exact_extract(loss_cropped, country_sf, "sum")
      country_biom <- exact_extract(biom_cropped, country_sf, "sum")
      
      # Convert to data frame
      country_data <- country_sf %>%
        st_drop_geometry() %>%
        mutate(
          loss_tc = country_loss,  # Already in tC (CORRECTED: removed *0.5)
          biomass_tc = country_biom
        ) %>%
        left_join(country_energy_shares, by = "country")
      
      # Calculate proportional loss for each country
      country_data <- country_data %>%
        mutate(
          prop_loss = if_else(biomass_tc > 0, loss_tc / biomass_tc, 0)
          #,prop_loss = pmin(prop_loss, 0.3)  # Cap at 30%
        )
      
      # Get total EU harvest volume for this year
      harvest_volume <- harvest_scenario %>%
        filter(year == !!year) %>%
        pull(harvest_m3)
      
      if (length(harvest_volume) == 0) next
      
      # Estimate country-level harvest based on biomass proportion
      total_biomass <- sum(country_data$biomass_tc, na.rm = TRUE)
      country_data <- country_data %>%
        mutate(
          harvest_share = biomass_tc / total_biomass,
          harvest_before = harvest_share * harvest_volume,
          harvest_after = harvest_before * (1 - prop_loss)
        )
      
      # Calculate energy biomass with country-specific shares
      country_data <- country_data %>%
        mutate(
          energy_before = harvest_before * conversion_factor * energy_share,
          energy_after = harvest_after * conversion_factor * energy_share
        )
      
      # Aggregate to EU27+UK level
      eu_energy_before <- sum(country_data$energy_before, na.rm = TRUE)
      eu_energy_after <- sum(country_data$energy_after, na.rm = TRUE)
      eu_reduction_pct <- if_else(eu_energy_before > 0, (1 - eu_energy_after/eu_energy_before) * 100, 0)
      
      # Store results
      results <- bind_rows(results, tibble(
        scenario = scenario,
        year = year,
        total_biomass_tc = total_biomass,
        fire_loss_tc = sum(country_data$loss_tc, na.rm = TRUE),
        harvest_before_m3 = harvest_volume,
        harvest_after_m3 = sum(country_data$harvest_after, na.rm = TRUE),
        energy_before_tc = eu_energy_before,
        energy_after_tc = eu_energy_after,
        reduction_pct = eu_reduction_pct,
        prop_loss = weighted.mean(country_data$prop_loss, country_data$biomass_tc, na.rm = TRUE)
      ))
    }, error = function(e) {
      message(glue("Error processing {scenario} {year}: {e$message}"))
    })
  }
  return(results)
}

# Process all scenarios
energy_data <- map_dfr(scenarios, calculate_fire_impact)

# 4. SANITY CHECKS ───────────────────────────────────────────────────────────
# 4.1 Check that fire loss doesn't exceed biomass
if ("fire_loss_tc" %in% names(energy_data) && "total_biomass_tc" %in% names(energy_data)) {
  biomass_check <- energy_data %>%
    filter(fire_loss_tc > total_biomass_tc)
  
  if (nrow(biomass_check) > 0) {
    message("CRITICAL ERROR: Fire loss exceeds total biomass in some years:")
    print(biomass_check)
    # Cap fire loss at biomass amount
    energy_data <- energy_data %>%
      mutate(fire_loss_tc = pmin(fire_loss_tc, total_biomass_tc))
  } else {
    message("✓ Fire loss never exceeds total biomass")
  }
} else {
  message("WARNING: Required columns missing for biomass check")
}

# 4.2 Check against expected ranges
range_checks <- energy_data %>%
  mutate(
    harvest_outside = harvest_before_m3 < expected_ranges$harvest_m3["min"] | 
      harvest_before_m3 > expected_ranges$harvest_m3["max"],
    biomass_outside = total_biomass_tc < expected_ranges$biomass_tc["min"] | 
      total_biomass_tc > expected_ranges$biomass_tc["max"],
    energy_outside = energy_before_tc < expected_ranges$energy_tc["min"] | 
      energy_before_tc > expected_ranges$energy_tc["max"],
    reduction_outside = reduction_pct < expected_ranges$reduction_pct["min"] | 
      reduction_pct > expected_ranges$reduction_pct["max"]
  )

if (any(range_checks$harvest_outside, na.rm = TRUE)) {
  message("WARNING: Harvest volumes outside expected range:")
  print(range_checks %>% filter(harvest_outside))
}

if (any(range_checks$biomass_outside, na.rm = TRUE)) {
  message("WARNING: Biomass values outside expected range:")
  print(range_checks %>% filter(biomass_outside))
}

if (any(range_checks$energy_outside, na.rm = TRUE)) {
  message("WARNING: Energy biomass outside expected range:")
  print(range_checks %>% filter(energy_outside))
}

if (any(range_checks$reduction_outside, na.rm = TRUE)) {
  message("WARNING: Reduction percentages outside expected range:")
  print(range_checks %>% filter(reduction_outside))
}

# 5. Save results ─────────────────────────────────────────────────────────────
write_csv(energy_data, file.path(results_dir, "Fire_Impact_Results.csv"))

# 6. Enhanced Creative Visualizations ───────────────────────────────────────
# 6.1 Harvest Volume Impact: Line Plot with Ribbon
p_harvest_line <- energy_data %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = harvest_before_m3 / 1e6, color = "Before Fire"), size = 1.2) +
  geom_line(aes(y = harvest_after_m3 / 1e6, color = "After Fire"), size = 1.2) +
  geom_ribbon(aes(ymin = harvest_after_m3 / 1e6, 
                  ymax = harvest_before_m3 / 1e6,
                  fill = "Fire Impact"), alpha = 0.3) +
  facet_wrap(~ scenario, ncol = 1,
             labeller = labeller(scenario = c(
               "Mit2p6" = "2.6°C Scenario",
               "Mit4p5" = "4.5°C Scenario",
               "Mit7p0" = "7.0°C Scenario"
             ))) +
  scale_color_manual(
    values = c("Before Fire" = "#1a9850", "After Fire" = "#d73027"),
    name = NULL
  ) +
  scale_fill_manual(values = c("Fire Impact" = "#fdae61"), name = NULL) +
  labs(
    title = "Harvest Volume Before and After Wildfires",
    subtitle = "EU27+UK projections (2015-2100)",
    x = "Year",
    y = "Harvest Volume (million m³)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# 6.2 Energy Biomass Impact: Streamgraph
p_energy_stream <- energy_data %>%
  select(year, scenario, energy_before_tc, energy_after_tc) %>%
  pivot_longer(cols = c(energy_before_tc, energy_after_tc), 
               names_to = "period", values_to = "energy") %>%
  mutate(period = factor(period, 
                         levels = c("energy_before_tc", "energy_after_tc"),
                         labels = c("Before Fire", "After Fire"))) %>%
  ggplot(aes(x = year, y = energy / 1e6, fill = period)) +
  geom_area(position = "stack", alpha = 0.85) +
  facet_wrap(~ scenario, ncol = 3,
             labeller = labeller(scenario = c(
               "Mit2p6" = "2.6°C Scenario",
               "Mit4p5" = "4.5°C Scenario",
               "Mit7p0" = "7.0°C Scenario"
             ))) +
  scale_fill_manual(
    values = c("Before Fire" = "#1a9850", "After Fire" = "#d73027"),
    name = "Energy Biomass"
  ) +
  labs(
    title = "Energy Biomass Availability Before and After Wildfires",
    subtitle = "EU27+UK projections across scenarios (2015-2100)",
    x = "Year",
    y = "Energy Biomass (MtC)",
    caption = "Using country-specific energy share factors"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# 6.3 Reduction Percentage: Animated Scatter Plot
# Identify key years with max reduction
key_years <- energy_data %>%
  group_by(scenario) %>%
  slice_max(reduction_pct, n = 1) %>%
  ungroup()

p_reduction_scatter <- ggplot(energy_data, aes(x = year, y = reduction_pct, color = scenario)) +
  geom_point(aes(size = reduction_pct), alpha = 0.7) +
  geom_line(aes(group = scenario), alpha = 0.5) +
  geom_label_repel(
    data = key_years,
    aes(label = paste0(round(reduction_pct, 1), "%")),
    size = 4,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  scale_color_viridis_d(
    option = "inferno",
    name = "Scenario",
    labels = c("Mit2p6" = "2.6°C", "Mit4p5" = "4.5°C", "Mit7p0" = "7.0°C")
  ) +
  scale_size_continuous(range = c(1, 8), name = "Reduction (%)") +
  labs(
    title = "Percentage Reduction in Energy Biomass Due to Wildfires",
    subtitle = "EU27+UK projections (2015-2100) with key years highlighted",
    x = "Year",
    y = "Reduction (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    panel.grid.minor = element_blank()
  )

# 6.4 Cumulative Impact: Radial Plot
cumulative_data <- energy_data %>%
  group_by(scenario) %>%
  summarise(
    cumulative_loss = sum(energy_before_tc - energy_after_tc, na.rm = TRUE) / 1e6,
    max_reduction = max(reduction_pct, na.rm = TRUE),
    years_affected = sum(reduction_pct > 0, na.rm = TRUE),
    .groups = "drop"
  )

p_cumulative_radial <- ggplot(cumulative_data, aes(x = scenario, y = cumulative_loss, fill = scenario)) +
  geom_col(width = 0.8, alpha = 0.85) +
  coord_polar() +
  geom_text(aes(label = paste0(round(cumulative_loss, 1), " MtC")), 
            position = position_stack(vjust = 0.5), size = 5, color = "white") +
  scale_fill_viridis_d(
    option = "viridis",
    name = "Scenario",
    labels = c("Mit2p6" = "2.6°C", "Mit4p5" = "4.5°C", "Mit7p0" = "7.0°C")
  ) +
  labs(
    title = "Cumulative Energy Biomass Loss Due to Wildfires",
    subtitle = "Total impact across scenarios (2015-2100)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# 6.5 NEW: Combined Harvest and Energy Biomass Volume Plot (2015-2070) ───────
# Compute energy biomass in volume units (Mm³)
energy_data <- energy_data %>%
  mutate(
    energy_before_Mm3 = energy_before_tc / (conversion_factor * 1e6),
    energy_after_Mm3 = energy_after_tc / (conversion_factor * 1e6)
  )

# Filter years to 2015-2070
energy_data_2015_2070 <- energy_data %>%
  filter(year <= 2070)

# Create the combined plot
p_combined_volume <- energy_data_2015_2070 %>%
  mutate(
    harvest_before_Mm3 = harvest_before_m3 / 1e6,
    harvest_after_Mm3 = harvest_after_m3 / 1e6
  ) %>%
  select(scenario, year, 
         harvest_before_Mm3, harvest_after_Mm3,
         energy_before_Mm3, energy_after_Mm3) %>%
  pivot_longer(
    cols = -c(scenario, year),
    names_to = "variable",
    values_to = "volume"
  ) %>%
  mutate(
    period = if_else(str_detect(variable, "before"), "Before Fire", "After Fire"),
    type = if_else(str_detect(variable, "harvest"), "Total Harvest", "Energy Biomass"),
    # Create a grouping variable for axis breaks
    axis_group = if_else(type == "Total Harvest", "High", "Low")
  ) %>%
  # Create a fake year offset for the low values to separate them visually
  mutate(
    year_offset = if_else(axis_group == "Low", year + 0.3, year)
  ) %>%
  ggplot(aes(x = year_offset, y = volume, color = type, group = interaction(type, period))) +
  
  # Shaded area to represent burnt biomass
  geom_ribbon(
    data = . %>% 
      group_by(scenario, year, type, axis_group) %>% 
      filter(n() == 2) %>% 
      arrange(desc(period)) %>% 
      reframe(
        x = year_offset,
        ymin = min(volume),
        ymax = max(volume)
      ),
    aes(x = x, ymin = ymin, ymax = ymax, fill = "Burnt Biomass"),
    alpha = 0.3, color = NA, inherit.aes = FALSE
  ) +
  
  # Main lines
  geom_line(aes(linetype = period), size = 1.8, alpha = 0.95) +
  
  # Facet with axis breaks illusion
  facet_grid(axis_group ~ scenario, scales = "free_y", space = "free_y",
             labeller = labeller(
               scenario = c(
                 "Mit2p6" = "Mitigation 2.6°C",
                 "Mit4p5" = "Mitigation 4.5°C",
                 "Mit7p0" = "Mitigation 7.0°C"
               ),
               axis_group = c("High" = "", "Low" = "")
             )) +
  
  # Scales and colors
  scale_color_manual(
    values = c("Total Harvest" = "#1a9850", "Energy Biomass" = "#4575b4"),
    name = "Volume Type"
  ) +
  scale_fill_manual(
    values = c("Burnt Biomass" = "#fc8d59"),
    name = "Fire Impact"
  ) +
  scale_linetype_manual(
    values = c("Before Fire" = "solid", "After Fire" = "11"),
    name = "Period"
  ) +
  scale_y_continuous(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_x_continuous(
    breaks = seq(2020, 2070, by = 10),
    limits = c(2015, 2070),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  # Titles and labels
  labs(
    title = "FIRE IMPACT ON FOREST BIOMASS RESOURCES (2015-2070)",
    subtitle = "EU27+UK projections of harvest volume and bioenergy feedstock availability",
    x = "Year",
    y = "Volume (Mm³/year)",
    caption = "Shaded areas represent biomass lost to wildfires"
  ) +
  
  # Theme with creative elements
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(1, "cm"),
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5, 
                              margin = margin(b = 10), color = "#d73027"),
    plot.subtitle = element_text(size = 18, hjust = 0.5, margin = margin(b = 20)),
    plot.caption = element_text(hjust = 0.5, size = 12, color = "gray40"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray80"),
    strip.text.x = element_text(face = "bold", size = 12, color = "black"),
    strip.text.y = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
  ) +
  
  # Custom legend with icons
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = 15, linetype = "solid"),
      title.position = "top",
      title.hjust = 0.5,
      order = 1
    ),
    linetype = guide_legend(
      override.aes = list(size = 1.5, color = "black"),
      title.position = "top",
      title.hjust = 0.5,
      order = 2
    ),
    fill = guide_legend(
      override.aes = list(alpha = 0.3, color = NA),
      title.position = "top",
      title.hjust = 0.5,
      order = 3
    )
  ) +
  
  # Decade markers
  geom_vline(
    xintercept = seq(2020, 2070, by = 10),
    color = "gray80",
    linetype = "dashed",
    alpha = 0.7
  )

# 7. Save Enhanced Visualizations ────────────────────────────────────────────
# Save plots
ggsave(file.path(results_dir, "Harvest_Volume_Impact.png"), 
       p_harvest_line, width = 12, height = 10, dpi = 300)

ggsave(file.path(results_dir, "Reduction_Percentage_Scatter.png"), 
       p_reduction_scatter, width = 12, height = 8, dpi = 300)
ggsave(file.path(results_dir, "Cumulative_Loss_Radial.png"), 
       p_cumulative_radial, width = 10, height = 10, dpi = 300)

# Save the new combined visualization
ggsave(file.path(results_dir, "Combined_Harvest_Energy_Volume_2015-2070.png"), 
       p_combined_volume, width = 18, height = 10, dpi = 300)

# Create combined report
combined_plot <- (p_harvest_line + p_energy_stream) / (p_reduction_scatter + p_cumulative_radial) +
  plot_annotation(
    title = "Comprehensive Fire Impact Analysis on Harvested Biomass",
    subtitle = "EU27+UK projections across climate scenarios (2015-2100)",
    caption = paste("Conversion factor:", conversion_factor, "tC/m³ | Country-specific energy shares"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5)
    )
  )

ggsave(file.path(results_dir, "Combined_Fire_Impact_Analysis.png"), 
       combined_plot, width = 18, height = 16, dpi = 300)

# Save data tables
energy_data %>%
  write_csv(file.path(results_dir, "Detailed_Fire_Impact_Results.csv"))

# Summary table
energy_data %>%
  group_by(scenario) %>%
  summarise(
    Total_Fire_Loss_MtC = sum(fire_loss_tc, na.rm = TRUE) / 1e6,
    Avg_Harvest_Reduction_Pct = mean(reduction_pct, na.rm = TRUE),
    Max_Harvest_Reduction_Pct = max(reduction_pct, na.rm = TRUE),
    Total_Energy_Loss_MtC = sum(energy_before_tc - energy_after_tc, na.rm = TRUE) / 1e6,
    Years_With_Significant_Fire = sum(prop_loss > 0.1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  write_csv(file.path(results_dir, "Summary_Fire_Impact_Results.csv"))

message("Analysis complete! Results saved to: ", results_dir)