# Enhanced European Wildfire Impacts Visualization
#   • Completely redesigned Plot 5 with creative visualization
#   • high resolution output
# ──────────────────────────────────────────────────────────────────────────────

library(terra)
library(tidyverse)
library(glue)
library(rnaturalearth)
library(sf)
library(viridis)
library(patchwork)
library(ggtext)
library(ggridges)
library(scales)
library(forcats)
library(ggrepel)

# Configuration
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
output_dir <- "H:/UPTAKE/Biomass_Loss_Results"
plot_dir <- file.path(output_dir, "Visualizations")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# European boundaries
europe <- ne_countries(scale = "medium", continent = "europe", returnclass = "sf") %>%
  st_crop(xmin = -15, xmax = 35, ymin = 35, ymax = 70)

# Palette setup
scenario_pal <- c("Mit2p6" = "#1a9850", "Mit4p5" = "#fee08b", "Mit7p0" = "#d73027")
forest_pal <- c("Afforested" = "#7fbc41", "Managed" = "#4d9221")

# 1. Load and prepare data with unit corrections ──────────────────────────────
load_annual_data <- function() {
  map_dfr(scenarios, function(sc) {
    read_csv(file.path(output_dir, sc, "annual_results.csv")) %>%  
      mutate(
        scenario = factor(sc, levels = c("Mit2p6", "Mit4p5", "Mit7p0")),
        # Apply unit corrections:
        across(c(total_burn, aff_burn, mgd_burn), ~ .x ),  
        across(c(total_loss, aff_loss, mgd_loss), ~ .x / 100)  
      )
  })
}

load_spatial_data <- function(year, variable) {
  map_dfr(scenarios, function(sc) {
    r <- rast(file.path(output_dir, sc, glue("{variable}_{year}.tif")))  
    
    # Apply unit corrections based on variable type
    if (str_detect(variable, "Burn")) {
      r <- r / 100  # Convert ha → km² for burnt area
    } else {
      r <- r / 100  # Remove erroneous *100 multiplier for loss
    }
    
    as.data.frame(r, xy = TRUE, na.rm = FALSE) %>%
      setNames(c("x", "y", "value")) %>%
      mutate(scenario = sc, year = year, variable = variable)
  })
}

# Load corrected annual results
annual_data <- load_annual_data()

# Load spatial data for key years
spatial_years <- c(2020, 2050, 2080, 2100)
spatial_data <- map_dfr(spatial_years, function(yr) {
  bind_rows(
    load_spatial_data(yr, "Total_Burn"),
    load_spatial_data(yr, "Total_Loss")
  )
})

# 2. Plot 1: Cumulative Loss with Shadow Area Under Lines ────────────────────
cumulative_loss <- annual_data %>%
  group_by(scenario) %>%
  mutate(cumul_loss = cumsum(total_loss)/1e6) %>%  # Convert to MtC
  ungroup()

p1 <- ggplot(cumulative_loss, aes(x = year, y = cumul_loss, color = scenario, fill = scenario)) +
  geom_area(position = "identity", alpha = 0.3) +
  geom_line(size = 1.5) +
  scale_color_manual(values = scenario_pal) +
  scale_fill_manual(values = scenario_pal) +
  scale_y_continuous(
    limits = c(0, max(cumulative_loss$cumul_loss) * 1.1),
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::comma
  ) +
  labs(title = "Cumulative Biomass Loss",
       subtitle = "Accumulated loss over time with shadow area (MtC)",
       x = "Year", y = "Cumulative Loss (MtC)",
       color = "Scenario", fill = "Scenario") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )

# 3. Plot 2: Enhanced Spatial-Temporal Evolution with Fixed NA Values ───────
spatial_plot_data <- spatial_data %>%
  filter(variable == "Total_Loss") %>%
  mutate(
    decade = (year %/% 10) * 10,
    value = value / 1e3,  # Convert to ktC (thousand tons of carbon)
    # Set NA values to 0 as requested
    value = ifelse(is.na(value), 0, value),
    # Create meaningful classes starting from >0
    value_class = cut(value, 
                      breaks = c(0, 0.001, 1, 5, 10, 50, 100, Inf),
                      labels = c("0", "0.001-1", "1-5", "5-10", "10-50", "50-100", "100+"),
                      include.lowest = TRUE)
  )

p2 <- ggplot(spatial_plot_data, aes(x = x, y = y)) +
  geom_tile(aes(fill = value_class)) +
  geom_sf(data = europe, fill = NA, color = "gray20", size = 0.3, inherit.aes = FALSE) +
  facet_grid(rows = vars(scenario), cols = vars(decade)) +
  scale_fill_viridis_d(
    option = "plasma",
    na.value = "transparent",
    name = "Biomass Loss (ktC)",
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +  # CORRECTED TO coord_sf
  labs(title = "Biomass Loss Evolution by Decade",
       subtitle = "Spatial patterns across scenarios with NA values set to 0") +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 10, hjust = 0.5)
  )

# 4. Plot 3: Burned Area vs. Loss Correlation ────────────────────────────────
p3_data <- annual_data %>%
  group_by(scenario, year) %>%
  summarise(
    total_burn = mean(total_burn),  # Now in km²
    total_loss = mean(total_loss)/1e6,  # Convert to MtC
    .groups = "drop"
  )

p3 <- ggplot(p3_data, aes(x = total_burn, y = total_loss, color = scenario)) +
  geom_point(aes(size = year), alpha = 0.7) +
  geom_path(aes(group = scenario), alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.7) +
  scale_color_manual(values = scenario_pal) +
  scale_size_continuous(range = c(1, 4), breaks = seq(2020, 2100, 20)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Burned Area vs. Biomass Loss Relationship",
       x = "Burned Area (km²)", 
       y = "Biomass Loss (MtC)",
       color = "Scenario",
       size = "Year") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right"
  )

# 5. Plot 4: Cumulative Impact Comparison ───────────────────────────────────
cumulative_summary <- annual_data %>%
  group_by(scenario) %>%
  summarise(
    total_burn = sum(total_burn, na.rm = TRUE),
    total_loss = sum(total_loss, na.rm = TRUE)/1e6,  # Convert to MtC
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(total_burn, total_loss), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "total_burn" = "Burned Area (km²)",
                         "total_loss" = "Biomass Loss (MtC)"))

p4 <- ggplot(cumulative_summary, aes(x = scenario, y = value, fill = scenario)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = scales::comma(round(value, 1))), 
            vjust = -0.5, size = 4, fontface = "bold") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = scenario_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Cumulative Impact (2015-2100)",
       subtitle = "Total burned area and biomass loss",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

# 6. NEW PLOT 5: Creative Forest Type Impact Comparison ──────────────────────
# Prepare data for creative visualization
forest_data <- annual_data %>%
  group_by(scenario) %>%
  summarise(
    aff_burn = sum(aff_burn, na.rm = TRUE),
    mgd_burn = sum(mgd_burn, na.rm = TRUE),
    aff_loss = sum(aff_loss, na.rm = TRUE)/1e6,  # Convert to MtC
    mgd_loss = sum(mgd_loss, na.rm = TRUE)/1e6,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(aff_burn, mgd_burn, aff_loss, mgd_loss),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    forest_type = ifelse(grepl("aff", metric), "Afforested", "Managed"),
    variable = ifelse(grepl("burn", metric), "Burned Area", "Biomass Loss"),
    unit = ifelse(grepl("burn", metric), "km²", "MtC"),
    # Create a combined label
    label = paste0(round(value, 1), " ", unit)
  )

# Create a completely new creative visualization - Bubble chart
# 6. Enhanced Streamgraph with Revised Color Scheme ────────────────────────
# Prepare data for streamgraph visualization
stream_data <- annual_data %>%
  select(year, scenario, aff_burn, mgd_burn, aff_loss, mgd_loss) %>%
  pivot_longer(
    cols = c(aff_burn, mgd_burn, aff_loss, mgd_loss),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    forest_type = ifelse(grepl("aff", metric), "Afforested", "Managed"),
    impact_type = ifelse(grepl("burn", metric), "Burned Area", "Biomass Loss"),
    # Create a combined variable for coloring
    impact_forest = paste(impact_type, forest_type, sep = " - "),
    # Normalize values for streamgraph
    value = ifelse(impact_type == "Biomass Loss", value / 1e6, value)  # Convert to MtC for loss
  ) %>%
  group_by(year, scenario, impact_type) %>%
  mutate(
    total_value = sum(value),
    percentage = value / total_value * 100
  ) %>%
  ungroup()

# Create a unified color palette with orange/red for Burned Area and green for Biomass Loss
impact_forest_pal <- c(
  "Burned Area - Afforested" = "#FFA07A",  
  "Burned Area - Managed" = "#FF4500",     
  "Biomass Loss - Afforested" = "seagreen3",
  "Biomass Loss - Managed" = "seagreen4"    
)

# Create a streamgraph visualization with impact-specific coloring
p5 <- ggplot(stream_data, aes(x = year, y = value, fill = impact_forest)) +
  geom_area(position = "fill", alpha = 0.85) +
  geom_line(aes(color = impact_forest), position = "fill", size = 0.5, alpha = 0.7) +
  facet_grid(impact_type ~ scenario, scales = "free_y") +
  scale_fill_manual(values = impact_forest_pal) +
  scale_color_manual(values = impact_forest_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Forest Type Impact Streamgraph",
    subtitle = "Temporal evolution of forest type contributions to impacts",
    x = "Year",
    y = "Percentage Contribution",
    fill = "Forest Type",
    color = "Forest Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1, "cm")
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 2, byrow = TRUE)
  )
# 7. Combine and save ULTRA-HIGH-QUALITY plots ──────────────────────────────
design <- "
AABB
AABB
CCDD
EEEE
EEEE
"

combined_plot <- p1 + p2 + p3 + p4 + p5 + 
  plot_layout(design = design) +
  plot_annotation(
    title = "Enhanced European Wildfire Impacts Analysis",
    subtitle = "2015-2100 Projections | FLAM-G4M Integration",
    caption = "Data: Processed biomass loss outputs | Visualization: Enhanced analysis",
    theme = theme(
      plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5, color = "gray30"),
      plot.caption = element_text(size = 14, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Save at maximum quality with increased size
ggsave(
  file.path(plot_dir, "Combined_Enhanced_Visualizations.png"),
  combined_plot,
  width = 24,  # Increased width
  height = 20, # Increased height
  dpi = 600,   # Ultra-high resolution
  bg = "white"
)

# 8. Save individual ultra-high-quality plots ────────────────────────────────
save_plot <- function(plot, name) {
  ggsave(
    file.path(plot_dir, name),
    plot,
    width = 14,  # Increased size
    height = 12, # Increased size
    dpi = 600,   # Ultra-high resolution
    bg = "white"
  )
}

save_plot(p1, "1_Cumulative_Loss_With_Shadow.png")
save_plot(p2, "2_Enhanced_Spatial_Evolution.png")
save_plot(p3, "3_Burn_Loss_Correlation.png")
save_plot(p4, "4_Cumulative_Impact.png")
save_plot(p5, "5_Creative_Forest_Type_Impact.png")

message("Enhanced visualizations saved to: ", plot_dir)








# CSV Extraction Script for Visualization Data
#   • Extracts data used in each visualization
#   • Saves as CSV files for reproducibility and analysis
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(terra)

# Configuration
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
output_dir <- "H:/UPTAKE/Biomass_Loss_Results"
csv_dir <- file.path(output_dir, "Visualization_Data_CSV")
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load annual data with unit corrections ──────────────────────────────────
load_annual_data <- function() {
  map_dfr(scenarios, function(sc) {
    read_csv(file.path(output_dir, sc, "annual_results.csv")) %>%  
      mutate(
        scenario = factor(sc, levels = c("Mit2p6", "Mit4p5", "Mit7p0")),
        # Apply unit corrections:
        across(c(total_burn, aff_burn, mgd_burn), ~ .x ),  
        across(c(total_loss, aff_loss, mgd_loss), ~ .x / 100)  
      )
  })
}

annual_data <- load_annual_data()

# 2. Extract Plot 1 Data: Cumulative Loss ───────────────────────────────────
cumulative_loss <- annual_data %>%
  group_by(scenario) %>%
  mutate(cumul_loss = cumsum(total_loss)/1e6) %>%  # Convert to MtC
  ungroup() %>%
  select(year, scenario, cumul_loss)

write_csv(cumulative_loss, file.path(csv_dir, "Plot1_Cumulative_Loss.csv"))

# 3. Extract Plot 2 Data: Spatial Evolution ──────────────────────────────────
extract_spatial_data <- function(year, variable) {
  map_dfr(scenarios, function(sc) {
    r <- rast(file.path(output_dir, sc, glue("{variable}_{year}.tif")))  
    
    # Apply unit corrections
    if (str_detect(variable, "Burn")) {
      r <- r / 100  # Convert ha → km²
    } else {
      r <- r / 100  # Remove erroneous *100 multiplier
    }
    
    as.data.frame(r, xy = TRUE, na.rm = FALSE) %>%
      setNames(c("x", "y", "value")) %>%
      mutate(scenario = sc, year = year, variable = variable)
  })
}

spatial_years <- c(2020, 2050, 2080, 2100)
spatial_data <- map_dfr(spatial_years, function(yr) {
  bind_rows(
    extract_spatial_data(yr, "Total_Burn"),
    extract_spatial_data(yr, "Total_Loss")
  )
}) %>%
  mutate(
    value = ifelse(variable == "Total_Loss", value / 1e3, value),  # Convert to ktC
    value = ifelse(is.na(value), 0, value)  # Set NA to 0
  )

write_csv(spatial_data, file.path(csv_dir, "Plot2_Spatial_Data.csv"))

# 4. Extract Plot 3 Data: Burned Area vs. Loss Correlation ───────────────────
correlation_data <- annual_data %>%
  group_by(scenario, year) %>%
  summarise(
    total_burn = mean(total_burn),  # km²
    total_loss = mean(total_loss)/1e6,  # MtC
    .groups = "drop"
  )

write_csv(correlation_data, file.path(csv_dir, "Plot3_Correlation_Data.csv"))

# 5. Extract Plot 4 Data: Cumulative Impact Comparison ─────────────────────
cumulative_impact <- annual_data %>%
  group_by(scenario) %>%
  summarise(
    total_burn = sum(total_burn, na.rm = TRUE),
    total_loss = sum(total_loss, na.rm = TRUE)/1e6,  # MtC
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(total_burn, total_loss), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "total_burn" = "Burned Area (km²)",
                         "total_loss" = "Biomass Loss (MtC)"))

write_csv(cumulative_impact, file.path(csv_dir, "Plot4_Cumulative_Impact.csv"))

# 6. Extract Plot 5 Data: Forest Type Impact Streamgraph ─────────────────────
streamgraph_data <- annual_data %>%
  select(year, scenario, aff_burn, mgd_burn, aff_loss, mgd_loss) %>%
  pivot_longer(
    cols = c(aff_burn, mgd_burn, aff_loss, mgd_loss),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    forest_type = ifelse(grepl("aff", metric), "Afforested", "Managed"),
    impact_type = ifelse(grepl("burn", metric), "Burned Area", "Biomass Loss"),
    value = ifelse(impact_type == "Biomass Loss", value / 1e6, value)  # Convert to MtC
  ) %>%
  group_by(year, scenario, impact_type) %>%
  mutate(
    total_value = sum(value),
    percentage = value / total_value * 100
  ) %>%
  ungroup() %>%
  select(year, scenario, forest_type, impact_type, value, percentage)

write_csv(streamgraph_data, file.path(csv_dir, "Plot5_Streamgraph_Data.csv"))

message("All visualization data saved to: ", csv_dir)
