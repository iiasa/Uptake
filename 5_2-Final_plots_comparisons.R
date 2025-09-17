# ──────────────────────────────────────────────────────────────────────────────
# Biomass Loss Analysis with UKESM (main) and GFDL (lower bound)
#   • Enhanced Visualizations with Consistent Color Schemes
#   • Optimized Plot Layouts
#   • Added Combined Annual + Cumulative Plot (Two Panels)
# ──────────────────────────────────────────────────────────────────────────────

library(terra)
library(tidyverse)
library(glue)
library(rnaturalearth)
library(sf)
library(viridis)
library(patchwork)
library(scales)
library(ggthemes)
library(ggtext)

# ──── CONFIGURATION ────────────────────────────────────────────────────────
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
years <- 2015:2100

# Directory setup
output_dir <- "H:/UPTAKE/Biomass_Loss_Results_UKESM"  # UKESM results (main)
gfdl_dir <- "H:/UPTAKE/Biomass_Loss_Results"          # GFDL results (lower bound)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create scenario mapping
scenario_map <- tibble(
  g4m_scenario = scenarios,
  flam_scenario = c("SSP126", "SSP245", "SSP370")
)

# Consistent color palette
palette_main <- c("Mit2p6" = "#1f77b4", "Mit4p5" = "#ff7f0e", "Mit7p0" = "#2ca02c")
palette_range <- c("Mit2p6" = "#aec7e8", "Mit4p5" = "#ffbb78", "Mit7p0" = "#98df8a")

# Font settings for high-quality output
font_size_base <- 16
font_size_title <- 20
font_size_subtitle <- 18

# ──── LOAD AND PROCESS DATA ────────────────────────────────────────────────
# Function to load results for a scenario
load_scenario_data <- function(scenario) {
  ukesm_file <- file.path(output_dir, scenario, "annual_results.csv")
  gfdl_file <- file.path(gfdl_dir, scenario, "annual_results.csv")
  
  if (!file.exists(ukesm_file) | !file.exists(gfdl_file)) {
    return(NULL)
  }
  
  ukesm <- read_csv(ukesm_file) %>%
    select(year, scenario, total_loss, total_burn) %>%
    rename(u_total_loss = total_loss, u_total_burn = total_burn)
  
  gfdl <- read_csv(gfdl_file) %>%
    select(year, total_loss) %>%
    rename(g_total_loss = total_loss)
  
  combined <- ukesm %>%
    left_join(gfdl, by = "year") %>%
    mutate(scenario = factor(scenario, levels = scenarios))
  
  return(combined)
}

# Load all data
all_data <- map_dfr(scenarios, load_scenario_data)

# Calculate cumulative loss
cumulative_data <- all_data %>%
  group_by(scenario) %>%
  arrange(year) %>%
  mutate(
    u_cumulative = cumsum(u_total_loss),
    g_cumulative = cumsum(g_total_loss)
  )

# Calculate relative change from baseline
baseline <- all_data %>%
  filter(year == 2015) %>%
  select(scenario, u_baseline = u_total_loss, g_baseline = g_total_loss)

relative_data <- all_data %>%
  left_join(baseline, by = "scenario") %>%
  mutate(
    u_relative = u_total_loss / u_baseline,
    g_relative = g_total_loss / g_baseline
  )

# ──── CREATE ENHANCED PLOTS ──────────────────────────────────────────────

# Custom theme for high-quality publications
theme_enhanced <- function(base_size = font_size_base) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = font_size_title, face = "bold", hjust = 0.5, 
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = font_size_subtitle, hjust = 0.5, 
                                   color = "gray30", margin = margin(b = 15)),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.position = "bottom",
      strip.text = element_text(size = base_size + 2, face = "bold", color = "black"),
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# 1. Temporal Trends with Model Range (Faceted by Scenario)
p1 <- ggplot(all_data, aes(x = year)) +
  geom_ribbon(aes(ymin = g_total_loss/1e8, ymax = u_total_loss/1e8, fill = scenario), 
              alpha = 0.3) +
  geom_line(aes(y = u_total_loss/1e8, color = scenario), size = 1.5) +
  geom_line(aes(y = g_total_loss/1e8), color = "gray40", linetype = "dashed", size = 1) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_color_manual(values = palette_main) +
  scale_fill_manual(values = palette_range) +
  labs(
    title = "Annual Biomass Loss Trends by Scenario",
    subtitle = "UKESM (solid) vs GFDL (dashed) with model range",
    x = "Year",
    y = "Annual Loss (MtC)",
    color = "Scenario",
    fill = "Model Range"
  ) +
  theme_enhanced() +
  theme(legend.position = "none")

# 2. Cumulative Loss Comparison (Faceted by Scenario)
p2 <- cumulative_data %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = g_cumulative/1e8, ymax = u_cumulative/1e8, fill = scenario), 
              alpha = 0.3) +
  geom_line(aes(y = u_cumulative/1e8, color = scenario), size = 1.5) +
  geom_line(aes(y = g_cumulative/1e8), color = "gray40", linetype = "dashed", size = 1) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_color_manual(values = palette_main) +
  scale_fill_manual(values = palette_range) +
  labs(
    title = "Cumulative Biomass Loss by Scenario",
    subtitle = "UKESM (solid) vs GFDL (dashed) with model range",
    x = "Year",
    y = "Cumulative Loss (MtC)"
  ) +
  theme_enhanced() +
  theme(legend.position = "none")

# 3. Combined Annual and Cumulative Plot (Two Panels)
# Left Panel: Annual Loss by Scenario
p_left <- ggplot(all_data, aes(x = year)) +
  geom_ribbon(aes(ymin = g_total_loss/1e8, ymax = u_total_loss/1e8, fill = scenario), 
              alpha = 0.3) +
  geom_line(aes(y = u_total_loss/1e8, color = scenario), size = 1.5) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_color_manual(values = palette_main) +
  scale_fill_manual(values = palette_range) +
  labs(
    title = "Annual Biomass Loss by Scenario",
    subtitle = "UKESM (solid) with GFDL range",
    x = "Year",
    y = "Annual Loss (MtC)"
  ) +
  theme_enhanced() +
  theme(legend.position = "none")

# Right Panel: Cumulative Loss All Scenarios Together
p_right <- cumulative_data %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = g_cumulative/1e8, ymax = u_cumulative/1e8, fill = scenario), 
              alpha = 0.3) +
  geom_line(aes(y = u_cumulative/1e8, color = scenario), size = 1.5) +
  scale_color_manual(values = palette_main) +
  scale_fill_manual(values = palette_range) +
  labs(
    title = "Cumulative Biomass Loss",
    subtitle = "All scenarios together with model ranges",
    x = "Year",
    y = "Cumulative Loss (MtC)"
  ) +
  theme_enhanced() +
  theme(legend.position = "bottom")

# Combine the two panels
p_combined <- p_left + p_right +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    title = "Biomass Loss Comparison: Annual vs Cumulative",
    subtitle = "Left: Annual Loss by Scenario | Right: Cumulative Loss All Scenarios",
    theme = theme(
      plot.title = element_text(size = font_size_title, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = font_size_subtitle, hjust = 0.5)
    )
  )

# 4. Relative Change from Baseline (Side-by-Side Comparison)
p3 <- ggplot(relative_data, aes(x = year)) +
  geom_ribbon(aes(ymin = g_relative, ymax = u_relative, fill = scenario), 
              alpha = 0.3) +
  geom_line(aes(y = u_relative, color = scenario), size = 1.5) +
  geom_line(aes(y = g_relative), color = "gray40", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray20", size = 1) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_color_manual(values = palette_main) +
  scale_fill_manual(values = palette_range) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    trans = "log10"
  ) +
  labs(
    title = "Relative Change in Biomass Loss from 2015 Baseline",
    subtitle = "UKESM (solid) vs GFDL (dashed) with model range",
    x = "Year",
    y = "Relative Change"
  ) +
  theme_enhanced() +
  theme(legend.position = "none")

# 5. Decadal Comparison of Average Loss (Faceted by Scenario)
decadal_data <- all_data %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(scenario, decade) %>%
  summarise(
    u_avg_loss = mean(u_total_loss, na.rm = TRUE),
    g_avg_loss = mean(g_total_loss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(decade >= 2020 & decade <= 2090)

p4 <- decadal_data %>%
  ggplot(aes(x = factor(decade), group = scenario)) +
  geom_col(aes(y = u_avg_loss/1e8, fill = scenario), alpha = 0.8) +
  geom_errorbar(aes(ymin = g_avg_loss/1e8, ymax = u_avg_loss/1e8),
                width = 0.3, color = "gray30", size = 1) +
  geom_point(aes(y = g_avg_loss/1e8), color = "gray20", size = 3) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_fill_manual(values = palette_main) +
  labs(
    title = "Average Decadal Biomass Loss by Scenario",
    subtitle = "UKESM (bars) vs GFDL (points) with model range",
    x = "Decade",
    y = "Average Annual Loss (MtC)"
  ) +
  theme_enhanced() +
  theme(legend.position = "none")

# 6. Model Agreement Analysis (Faceted by Scenario)
p5 <- all_data %>%
  mutate(model_diff = (u_total_loss - g_total_loss)/1e8) %>%
  ggplot(aes(x = year, y = model_diff, fill = scenario)) +
  geom_area(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_fill_manual(values = palette_main) +
  labs(
    title = "Model Difference Analysis (UKESM - GFDL)",
    subtitle = "Magnitude of difference between models by scenario",
    x = "Year",
    y = "Difference in Annual Loss (MtC)"
  ) +
  theme_enhanced() +
  theme(legend.position = "none")

# ──── SAVE OUTPUT ────────────────────────────────────────────
# Create comparison directory
comparison_dir <- file.path(output_dir, "Model_Comparison")
dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

# Save plots as PNG
ggsave(file.path(comparison_dir, "01_Temporal_Trends.png"), p1, 
       width = 10, height = 12, dpi = 300)
ggsave(file.path(comparison_dir, "02_Cumulative_Loss.png"), p2, 
       width = 10, height = 12, dpi = 300)
ggsave(file.path(comparison_dir, "03_Combined_Annual_Cumulative.png"), p_combined, 
       width = 16, height = 8, dpi = 300)
ggsave(file.path(comparison_dir, "04_Relative_Change.png"), p3, 
       width = 14, height = 8, dpi = 300)
ggsave(file.path(comparison_dir, "05_Decadal_Comparison.png"), p4, 
       width = 14, height = 8, dpi = 300)
ggsave(file.path(comparison_dir, "06_Model_Difference.png"), p5, 
       width = 14, height = 8, dpi = 300)

# Create combined report (PNG only)
combined_plot <- (p1 + plot_layout(tag_level = 'new')) / 
  (p2 + plot_layout(tag_level = 'new')) / 
  (p_combined + plot_layout(tag_level = 'new')) / 
  (p3 + plot_layout(tag_level = 'new')) / 
  (p4 + plot_layout(tag_level = 'new')) / 
  (p5 + plot_layout(tag_level = 'new')) +
  plot_annotation(
    title = "Comprehensive Biomass Loss Model Comparison",
    subtitle = "UKESM vs GFDL Analysis Across Scenarios",
    caption = "Data: UKESM (primary model) | GFDL (lower bound range)",
    theme = theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 22, hjust = 0.5),
      plot.caption = element_text(size = 16, hjust = 1)
    )
  )

ggsave(file.path(comparison_dir, "Full_Model_Comparison_Report.png"), 
       combined_plot, width = 16, height = 28, dpi = 300)

message("✅ All comparison plots created successfully!")
message("Results saved to: ", comparison_dir)