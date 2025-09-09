# Load required libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(dplyr)

# Read the CSV file
file_path <- "H:/UPTAKE/tabs_gui_ACCREU_v01_18062025_scen_csv_alln1.csv"
data <- fread(file_path,header=T)



# List of EU27 countries + UK
eu27_uk_countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark",
                       "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
                       "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
                       "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "United Kingdom")

# Filter for EU27 + UK countries and the relevant parameters
filtered_data <- data %>%
  filter(Country %in% eu27_uk_countries,
         Parameter %in% c("harvest_total_m3year", "residue_extract_m3"))

# Pivot to wide format for easier calculation
wide_data <- filtered_data %>%
  pivot_longer(cols = -c(Country, Scenario, Parameter),
               names_to = "Year",
               values_to = "Value") %>%
  mutate(Year = as.numeric(Year),
         Value = as.numeric(Value)) %>%
  pivot_wider(names_from = Parameter, values_from = Value)

# Calculate EU27+UK totals by Year and Scenario
eu_totals <- wide_data %>%
  group_by(Scenario, Year) %>%
  summarise(
    Total_Harvest = sum(harvest_total_m3year, na.rm = TRUE),
    Total_Residues = sum(residue_extract_m3, na.rm = TRUE),
    .groups = "drop"
  )

# Create the plot
ggplot(eu_totals, aes(x = Year)) +
  geom_line(aes(y = Total_Harvest, color = "Total Harvest"), size = 1) +
  geom_line(aes(y = Total_Residues, color = "Total Residues"), size = 1) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    title = "EU27+UK Total Harvest and Residue Extraction by Scenario",
    y = "Volume (m³)",
    x = "Year",
    color = "Parameter"
  ) +
  scale_color_manual(values = c("Total Harvest" = "blue", "Total Residues" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot
ggsave("EU_harvest_residues_analysis.png", width = 12, height = 8, dpi = 300)

# Display summary statistics
summary_stats <- eu_totals %>%
  group_by(Scenario) %>%
  summarise(
    Avg_Harvest = mean(Total_Harvest, na.rm = TRUE),
    Avg_Residues = mean(Total_Residues, na.rm = TRUE),
    Max_Harvest = max(Total_Harvest, na.rm = TRUE),
    Max_Residues = max(Total_Residues, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Write results to CSV
write_csv(eu_totals, "H:/UPTAKE/EU_harvest_residues_totals.csv")
write_csv(summary_stats, "H:/UPTAKE/EU_harvest_residues_summary.csv")

# Filter for specific scenarios
selected_scenarios <- c("Ref_Mit2p6__RCPref__NONE", "Ref_Mit4p5__RCPref__NONE", "Ref_Mit7p0__RCPref__NONE")

filtered_eu_totals <- eu_totals %>%
  filter(Scenario %in% selected_scenarios)

# Write filtered results to CSV
write_csv(filtered_eu_totals, "H:/UPTAKE/EU_harvest_residues_selected_scenarios.csv")

# Create a plot for the selected scenarios
ggplot(filtered_eu_totals, aes(x = Year)) +
  geom_line(aes(y = Total_Harvest, color = "Total Harvest"), size = 1) +
  geom_line(aes(y = Total_Residues, color = "Total Residues"), size = 1) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    title = "EU27+UK Total Harvest and Residue Extraction (Selected Scenarios)",
    y = "Volume (m³)",
    x = "Year",
    color = "Parameter"
  ) +
  scale_color_manual(values = c("Total Harvest" = "blue", "Total Residues" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Filter for specific scenarios
selected_scenarios <- c("Ref_Mit2p6__RCPref__NONE", "Ref_Mit4p5__RCPref__NONE", "Ref_Mit7p0__RCPref__NONE")
filtered_eu_totals <- eu_totals %>% filter(Scenario %in% selected_scenarios)


# Filter for specific scenarios
selected_scenarios <- c("Ref_Mit2p6__RCPref__NONE", "Ref_Mit4p5__RCPref__NONE", "Ref_Mit7p0__RCPref__NONE")
filtered_eu_totals <- eu_totals %>% filter(Scenario %in% selected_scenarios)

# Convert back to long format for easier plotting
plot_data <- filtered_eu_totals %>%
  pivot_longer(cols = c(Total_Harvest, Total_Residues),
               names_to = "Parameter",
               values_to = "Value")

# Preprocess scenario names and create facets
plot_data <- plot_data %>%
  mutate(
    Scenario = case_when(
      Scenario == "Ref_Mit2p6__RCPref__NONE" ~ "Mit2.6",
      Scenario == "Ref_Mit4p5__RCPref__NONE" ~ "Mit4.5",
      Scenario == "Ref_Mit7p0__RCPref__NONE" ~ "Mit7.0",
      TRUE ~ Scenario
    )
  )

# Create the faceted plot with shaded areas
ggplot(plot_data, aes(x = Year, y = Value, 
                      fill = Parameter, 
                      color = Parameter)) +
  # Add shaded areas first (underneath lines)
  geom_ribbon(aes(ymin = 0, ymax = Value), alpha = 0.5) +
  # Add lines on top of shaded areas
  geom_line(size = 1.2) +
  # Create separate panels for each scenario
  facet_wrap(~ Scenario, ncol = 3) +
  # Custom fill colors for ribbons (light tones)
  scale_fill_manual(
    values = c(
      "Total_Harvest" = "#E5F5E0",  # Light green for harvest
      "Total_Residues" = "#A1D99B"   # Light green for residues
    ),
    labels = c("Total Harvest", "Total Residues")
  ) +
  # Custom color palette for lines (dark tones)
  scale_color_manual(
    values = c(
      "Total_Harvest" = "#006D2C",  # Dark green for harvest
      "Total_Residues" = "#238B45"   # Dark green for residues
    ),
    labels = c("Total Harvest", "Total Residues")
  ) +
  # Labels and titles
  labs(
    title = "EU27+UK Forest Biomass Comparison Across Mitigation Scenarios",
    subtitle = "Comparing total harvest and residue extraction volumes",
    y = "Volume (million m³)",
    x = "Year",
    fill = "Biomass Type",
    color = "Biomass Type"
  ) +
  # Format y-axis in millions
  scale_y_continuous(labels = function(x) x / 1e6) +
  # Theme settings
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    panel.grid.minor = element_line(color = "grey90"),
    panel.grid.major = element_line(color = "grey85"),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 12, face = "bold")
  ) 
  # Direct labels at end of lines
 

# Save the plot
ggsave("H:/UPTAKE/EU_biomass_comparison_creative.png", width = 14, height = 8, dpi = 300)



# Create a comprehensive visualization with multiple elements
library(ggtext)
library(patchwork)

# Prepare data for plotting
plot_data <- eu_totals %>%
  filter(Scenario %in% selected_scenarios) %>%
  mutate(Scenario_clean = case_when(
    Scenario == "Ref_Mit2p6__RCPref__NONE" ~ "2.6°C Mitigation",
    Scenario == "Ref_Mit4p5__RCPref__NONE" ~ "4.5°C Mitigation",
    Scenario == "Ref_Mit7p0__RCPref__NONE" ~ "7.0°C Mitigation"
  )) %>%
  pivot_longer(cols = c(Total_Harvest, Total_Residues),
               names_to = "Biomass_Type",
               values_to = "Volume") %>%
  mutate(Biomass_Type = ifelse(Biomass_Type == "Total_Harvest", "Harvest", "Residues"),
         Volume_Millions = Volume / 1e6)

# Create main line plot
main_plot <- ggplot(plot_data, aes(x = Year, y = Volume_Millions,
                                   color = Scenario_clean,
                                   linetype = Biomass_Type)) +
  geom_line(size = 1.2, alpha = 0.9) +
  scale_color_manual(values = c("#1a9850", "#fd8d3c", "#d73027")) + # Green to red gradient
  scale_linetype_manual(values = c("solid", "longdash")) +
  labs(x = "Year", y = "Volume (million m³)",
       color = "Climate Scenario", linetype = "Biomass Type") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))

# Create small multiples for each scenario
scenario_plots <- plot_data %>%
  group_by(Scenario_clean) %>%
  nest() %>%
  mutate(plot = map2(data, Scenario_clean, ~{
    ggplot(.x, aes(x = Year, y = Volume_Millions, fill = Biomass_Type)) +
      geom_area(position = "identity", alpha = 0.7) +
      scale_fill_manual(values = c("#1a9850", "#a6d96a")) +
      labs(title = .y, x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 10, hjust = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }))

# Combine plots using patchwork
combined_plot <- main_plot +
  wrap_plots(scenario_plots$plot, nrow = 1) +
  plot_layout(ncol = 1, heights = c(3, 1)) +
  plot_annotation(
    title = "EU27+UK Forest Biomass Projections Under Different Climate Scenarios",
    subtitle = "Main plot shows trends across scenarios, while small multiples highlight composition differences",
    caption = "Data source: UPTAKE project | Harvest = Primary wood products | Residues = Secondary biomass",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
      plot.caption = element_text(size = 9, color = "grey50", hjust = 1)
    )
  )

# Add some annotations
final_plot <- combined_plot +
  geom_text(
    data = plot_data %>% filter(Year == max(Year)),
    aes(label = paste0(round(Volume_Millions, 0), "M"),
        x = Year + 2, y = Volume_Millions),
    size = 3, hjust = 0, color = "grey30"
  )

# Save the final plot
ggsave("EU_forest_biomass_creative_visualization.png",
       plot = final_plot,
       width = 14, height = 10, dpi = 300)

# Print the plot
final_plot