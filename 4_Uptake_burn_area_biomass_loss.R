# ──────────────────────────────────────────────────────────────────────────────
# Biomass Loss Analysis for Multiple Scenarios and Years
#   • Processes all years (2015-2100) and scenarios (Mit2p6, Mit4p5, Mit7p0)
#   • Calculates burnt area and biomass loss for afforested and managed forests
#   • Generates creative visualizations for temporal trends and spatial patterns
# ──────────────────────────────────────────────────────────────────────────────

library(terra)
library(tidyverse)
library(glue)
library(rnaturalearth)
library(sf)
library(viridis)
library(patchwork)
library(gganimate)
library(transformr)
library(scales)

# ──── CONFIGURATION ────────────────────────────────────────────────────────
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
flam_scenarios <- c("SSP126", "SSP245", "SSP370")
years <- 2015:2100

# Directory setup
flam_dir <- "H:/UPTAKE/FLAM_Output"
g4m_dir <- "H:/Uptake"
prop_dir <- "H:/UPTAKE/Area_Proportions"
output_dir <- "H:/UPTAKE/Biomass_Loss_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# European boundaries
europe <- ne_countries(scale = "medium", continent = "europe", returnclass = "sf") %>%
  st_crop(xmin = -15, xmax = 35, ymin = 35, ymax = 70)

# Create scenario mapping
scenario_map <- tibble(
  g4m_scenario = scenarios,
  flam_scenario = flam_scenarios
)

# ──── MAIN PROCESSING FUNCTION ────────────────────────────────────────────
process_scenario <- function(g4m_scenario, flam_scenario) {
  message("\n", strrep("═", 80))
  message("STARTING PROCESSING FOR SCENARIO: ", g4m_scenario)
  message(strrep("═", 80))

  scenario_out_dir <- file.path(output_dir, g4m_scenario)
  dir.create(scenario_out_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize results storage
  annual_results <- tibble()
  spatial_data <- list()

  # ──── PROCESS EACH YEAR ────────────────────────────────────────────────
  for (year in years) {
    message("\nProcessing year: ", year)

    # 1. FLAM BURNT AREA PROCESSING ───────────────────────────────────────
    flam_files <- list.files(
      file.path(flam_dir, glue("{flam_scenario}_A_burn_monthly")),
      pattern = glue("{year}.*\\.tif$"),
      full.names = TRUE
    )

    if (length(flam_files) == 0) {
      message("No FLAM files found for ", year)
      next
    }

    # Calculate annual burnt area
    burn_stack <- rast(flam_files)
    annual_burn <- sum(burn_stack)
    flam_ext <- ext(annual_burn)

    # 2. FOREST PROPORTIONS ───────────────────────────────────────────────
    aff_prop <- rast(file.path(prop_dir, g4m_scenario,
                               glue("{g4m_scenario}_afforested_proportion_stack.tif")))
    mgd_prop <- rast(file.path(prop_dir, g4m_scenario,
                               glue("{g4m_scenario}_managed_proportion_stack.tif")))

    aff_prop_yr <- aff_prop[[glue("y{year}")]]
    mgd_prop_yr <- mgd_prop[[glue("y{year}")]]

    # Crop to FLAM extent
    aff_prop_cropped <- crop(aff_prop_yr, flam_ext)
    mgd_prop_cropped <- crop(mgd_prop_yr, flam_ext)

    # 3. BURNT AREA BY FOREST TYPE ────────────────────────────────────────
    aff_burn <- annual_burn * aff_prop_cropped
    mgd_burn <- annual_burn * mgd_prop_cropped

    # Save burnt area rasters
    writeRaster(
      aff_burn,
      file.path(scenario_out_dir, glue("Afforested_Burn_{year}.tif")),
      overwrite = TRUE
    )
    writeRaster(
      mgd_burn,
      file.path(scenario_out_dir, glue("Managed_Burn_{year}.tif")),
      overwrite = TRUE
    )
    writeRaster(
      annual_burn,
      file.path(scenario_out_dir, glue("Total_Burn_{year}.tif")),
      overwrite = TRUE
    )

    # 4. BIOMASS DENSITY PROCESSING ───────────────────────────────────────
    aff_dens <- rast(file.path(g4m_dir, g4m_scenario, glue("{g4m_scenario}_Biom_af_{year}.tif")))
    mgd_dens <- rast(file.path(g4m_dir, g4m_scenario, glue("{g4m_scenario}_Biom_fm_total_{year}.tif")))

    # Convert to tC/ha
    cell_area_km2 <- cellSize(aff_dens, unit = "km")
    cell_area_ha <- cell_area_km2 * 100
    aff_dens_tcha <- aff_dens / cell_area_ha
    mgd_dens_tcha <- mgd_dens / cell_area_ha

    # Crop and resample
    aff_dens_crop <- crop(aff_dens_tcha, flam_ext)
    mgd_dens_crop <- crop(mgd_dens_tcha, flam_ext)
    aff_dens_res <- resample(aff_dens_crop, aff_burn, method = "bilinear")
    mgd_dens_res <- resample(mgd_dens_crop, mgd_burn, method = "bilinear")

    # 5. BIOMASS LOSS CALCULATION ────────────────────────────────────────
    aff_loss <- aff_burn * aff_dens_res * 100
    mgd_loss <- mgd_burn * mgd_dens_res * 100
    total_loss <- aff_loss + mgd_loss

    # Save loss rasters
    writeRaster(
      aff_loss,
      file.path(scenario_out_dir, glue("Afforested_Loss_{year}.tif")),
      overwrite = TRUE
    )
    writeRaster(
      mgd_loss,
      file.path(scenario_out_dir, glue("Managed_Loss_{year}.tif")),
      overwrite = TRUE
    )
    writeRaster(
      total_loss,
      file.path(scenario_out_dir, glue("Total_Loss_{year}.tif")),
      overwrite = TRUE
    )

    # 6. STORE RESULTS ───────────────────────────────────────────────────
    annual_results <- rbind(annual_results, tibble(
      year = year,
      scenario = g4m_scenario,
      total_burn = sum(values(annual_burn), na.rm = TRUE),
      aff_burn = sum(values(aff_burn), na.rm = TRUE),
      mgd_burn = sum(values(mgd_burn), na.rm = TRUE),
      aff_loss = sum(values(aff_loss), na.rm = TRUE),
      mgd_loss = sum(values(mgd_loss), na.rm = TRUE),
      total_loss = sum(values(total_loss), na.rm = TRUE)
    ))

    # Store spatial data for key years
    if (year %in% c(2015, 2030, 2050, 2070, 2100)) {
      spatial_data[[as.character(year)]] <- list(
        aff_loss = aff_loss,
        mgd_loss = mgd_loss,
        total_loss = total_loss
      )
    }

    message("✓ Completed year ", year)
  }

  # Save annual results
  write_csv(annual_results, file.path(scenario_out_dir, "annual_results.csv"))

  # ──── CREATE VISUALIZATIONS ─────────────────────────────────────────────
  message("\nCreating visualizations for ", g4m_scenario)

  # 1. TEMPORAL TRENDS PLOT ───────────────────────────────────────────────
  p_trend <- ggplot(annual_results, aes(x = year)) +
    geom_line(aes(y = total_loss/1e9, color = "Total Loss"), size = 1.2) +
    geom_line(aes(y = aff_loss/1e9, color = "Afforested"), size = 1) +
    geom_line(aes(y = mgd_loss/1e9, color = "Managed"), size = 1) +
    geom_area(aes(y = total_burn/max(total_burn)*max(total_loss/1e9)*0.8,
                  fill = "Burnt Area"), alpha = 0.3) +
    scale_color_manual(
      name = "Biomass Loss",
      values = c("Total Loss" = "#e41a1c", "Afforested" = "#4daf4a", "Managed" = "#377eb8")
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("Burnt Area" = "#ff7f00")
    ) +
    scale_y_continuous(
      name = "Biomass Loss (GtC)",
      sec.axis = sec_axis(~ . / max(annual_results$total_loss/1e9) * max(annual_results$total_burn),
                          name = "Burnt Area (km²)")
    ) +
    labs(
      title = glue("Wildfire Biomass Loss: {g4m_scenario}"),
      subtitle = "Annual trends with burnt area comparison",
      x = "Year"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    file.path(scenario_out_dir, "01_temporal_trends.png"),
    p_trend,
    width = 12,
    height = 8,
    dpi = 300
  )

  # 2. SPATIAL EVOLUTION PLOT ─────────────────────────────────────────────
  spatial_plots <- list()
  for (yr in names(spatial_data)) {
    # Create ocean mask
    europe_vect <- vect(europe)
    ocean_mask <- rasterize(europe_vect, spatial_data[[yr]]$total_loss, background = NA)

    # Apply ocean mask
    total_loss_plot <- mask(spatial_data[[yr]]$total_loss, ocean_mask, inverse = FALSE)

    # Convert to data frame
    loss_df <- as.data.frame(total_loss_plot, xy = TRUE, na.rm = FALSE) %>%
      setNames(c("x", "y", "loss"))

    # Create plot
    p_spatial <- ggplot() +
      geom_tile(data = loss_df, aes(x = x, y = y, fill = loss)) +
      geom_sf(data = europe, fill = NA, color = "gray30", size = 0.2, inherit.aes = FALSE) +
      scale_fill_viridis_c(
        option = "inferno",
        trans = "sqrt",
        name = "Biomass Loss (tC)",
        labels = scales::comma,
        na.value = "transparent"
      ) +
      coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
      labs(
        title = glue("Total Biomass Loss: {g4m_scenario} {yr}"),
        subtitle = "Spatial Distribution"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "bottom"
      )

    spatial_plots[[yr]] <- p_spatial
    ggsave(
      file.path(scenario_out_dir, glue("spatial_{yr}.png")),
      p_spatial,
      width = 10,
      height = 8,
      dpi = 300
    )
  }

  # Combine spatial plots
  if (length(spatial_plots) > 0) {
    p_spatial_grid <- wrap_plots(spatial_plots, ncol = 2) +
      plot_annotation(
        title = glue("Evolution of Biomass Loss: {g4m_scenario}"),
        subtitle = "Key years: 2015, 2030, 2050, 2070, 2100",
        theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
      )

    ggsave(
      file.path(scenario_out_dir, "02_spatial_evolution.png"),
      p_spatial_grid,
      width = 16,
      height = 12,
      dpi = 300
    )
  }

  # 3. ANIMATED SPATIAL EVOLUTION ─────────────────────────────────────────
  # Prepare data for animation
  # anim_data <- tibble()
  # for (yr in years) {
  #   loss_file <- file.path(scenario_out_dir, glue("Total_Loss_{yr}.tif"))
  #   if (!file.exists(loss_file)) next
  # 
  #   total_loss <- rast(loss_file)
  #   ocean_mask <- rasterize(europe_vect, total_loss, background = NA)
  #   total_loss_plot <- mask(total_loss, ocean_mask, inverse = FALSE)
  # 
  #   loss_df <- as.data.frame(total_loss_plot, xy = TRUE, na.rm = FALSE) %>%
  #     setNames(c("x", "y", "loss")) %>%
  #     mutate(year = yr, scenario = g4m_scenario)
  # 
  #   anim_data <- bind_rows(anim_data, loss_df)
  # }
  # 
  # if (nrow(anim_data) > 0) {
  #   p_anim <- ggplot() +
  #     geom_tile(data = anim_data, aes(x = x, y = y, fill = loss)) +
  #     geom_sf(data = europe, fill = NA, color = "gray30", size = 0.2, inherit.aes = FALSE) +
  #     scale_fill_viridis_c(
  #       option = "inferno",
  #       trans = "sqrt",
  #       name = "Biomass Loss (tC)",
  #       labels = scales::comma,
  #       na.value = "transparent"
  #     ) +
  #     coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
  #     labs(
  #       title = "Wildfire Biomass Loss: {closest_state}",
  #       subtitle = glue("Scenario: {g4m_scenario}"),
  #       caption = "Data: FLAM-G4M Integration"
  #     ) +
  #     theme_void() +
  #     theme(
  #       plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
  #       plot.subtitle = element_text(hjust = 0.5, size = 16)
  #     ) +
  #     transition_states(year, transition_length = 2, state_length = 1) +
  #     enter_fade() +
  #     exit_fade()
  # 
  #   animate(
  #     p_anim,
  #     nframes = 100,
  #     fps = 10,
  #     width = 1000,
  #     height = 800,
  #     renderer = gifski_renderer(file.path(scenario_out_dir, "03_spatial_evolution.gif"))
  #   )
  # }

  # 4. FOREST TYPE COMPARISON ─────────────────────────────────────────────
  p_forest <- annual_results %>%
    pivot_longer(cols = c(aff_loss, mgd_loss),
                 names_to = "forest_type", values_to = "loss") %>%
    mutate(forest_type = factor(forest_type,
                                levels = c("aff_loss", "mgd_loss"),
                                labels = c("Afforested", "Managed"))) %>%
    ggplot(aes(x = year, y = loss/1e9, fill = forest_type)) +
    geom_area(position = "stack", alpha = 0.8) +
    geom_line(aes(y = total_burn/max(total_burn)*max(total_loss/1e9)*0.8),
              color = "#ff7f00", size = 1.2) +
    scale_fill_manual(
      values = c("Afforested" = "#4daf4a", "Managed" = "#377eb8")
    ) +
    labs(
      title = glue("Biomass Loss by Forest Type: {g4m_scenario}"),
      subtitle = "Stacked area chart with burnt area trend",
      x = "Year",
      y = "Biomass Loss (GtC)",
      fill = "Forest Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    file.path(scenario_out_dir, "04_forest_type_comparison.png"),
    p_forest,
    width = 12,
    height = 8,
    dpi = 300
  )

  message("✓ Completed visualizations for ", g4m_scenario)
  return(annual_results)
}

# ──── PROCESS ALL SCENARIOS ───────────────────────────────────────────────
all_results <- map2_dfr(scenario_map$g4m_scenario,
                        scenario_map$flam_scenario,
                        process_scenario)

# ──── CREATE COMPARISON PLOTS ─────────────────────────────────────────────
message("\nCreating comparison plots across scenarios")

# 1. SCENARIO COMPARISON PLOT ──────────────────────────────────────────────
p_scenario <- all_results %>%
  group_by(scenario, year) %>%
  summarise(total_loss = sum(total_loss, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = total_loss/1e9, color = scenario)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Mit2p6" = "#1b9e77", "Mit4p5" = "#d95f02", "Mit7p0" = "#7570b3"),
    labels = c("Mit2p6" = "2.6°C", "Mit4p5" = "4.5°C", "Mit7p0" = "7.0°C")
  ) +
  labs(
    title = "Total Biomass Loss Across Climate Scenarios",
    subtitle = "Comparison of different mitigation pathways",
    x = "Year",
    y = "Total Biomass Loss (GtC)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "scenario_comparison.png"),
  p_scenario,
  width = 14,
  height = 9,
  dpi = 300
)

# 2. CUMULATIVE LOSS COMPARISON ────────────────────────────────────────────
p_cumulative <- all_results %>%
  group_by(scenario, year) %>%
  summarise(total_loss = sum(total_loss, na.rm = TRUE), .groups = "drop") %>%
  group_by(scenario) %>%
  mutate(cumulative_loss = cumsum(total_loss)/1e9) %>%
  ggplot(aes(x = year, y = cumulative_loss, color = scenario)) +
  geom_line(size = 1.5) +
  scale_color_manual(
    values = c("Mit2p6" = "#1b9e77", "Mit4p5" = "#d95f02", "Mit7p0" = "#7570b3"),
    labels = c("Mit2p6" = "2.6°C", "Mit4p5" = "4.5°C", "Mit7p0" = "7.0°C")
  ) +
  labs(
    title = "Cumulative Biomass Loss (2015-2100)",
    subtitle = "Total carbon lost to wildfires under different scenarios",
    x = "Year",
    y = "Cumulative Loss (GtC)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "cumulative_loss.png"),
  p_cumulative,
  width = 14,
  height = 9,
  dpi = 300
)

# 3. FINAL YEAR COMPARISON MAP ─────────────────────────────────────────────
final_year_data <- tibble()
for (scenario in scenarios) {
  loss_file <- file.path(output_dir, scenario, glue("Total_Loss_2100.tif"))
  if (!file.exists(loss_file)) next

  total_loss <- rast(loss_file)
  ocean_mask <- rasterize(vect(europe), total_loss, background = NA)
  total_loss_plot <- mask(total_loss, ocean_mask, inverse = FALSE)

  loss_df <- as.data.frame(total_loss_plot, xy = TRUE, na.rm = FALSE) %>%
    setNames(c("x", "y", "loss")) %>%
    mutate(scenario = scenario)

  final_year_data <- bind_rows(final_year_data, loss_df)
}

if (nrow(final_year_data) > 0) {
  p_final_year <- ggplot() +
    geom_tile(data = final_year_data, aes(x = x, y = y, fill = loss)) +
    geom_sf(data = europe, fill = NA, color = "gray30", size = 0.2, inherit.aes = FALSE) +
    facet_wrap(~ scenario, nrow = 2,
               labeller = labeller(scenario = c(
                 "Mit2p6" = "2.6°C Scenario",
                 "Mit4p5" = "4.5°C Scenario",
                 "Mit7p0" = "7.0°C Scenario"
               ))) +
    scale_fill_viridis_c(
      option = "inferno",
      trans = "sqrt",
      name = "Biomass Loss (tC)",
      labels = scales::comma,
      na.value = "transparent"
    ) +
    coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
    labs(
      title = "Projected Biomass Loss in 2100",
      subtitle = "Comparison across climate scenarios"
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 16),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  ggsave(
    file.path(output_dir, "final_year_comparison.png"),
    p_final_year,
    width = 16,
    height = 12,
    dpi = 300
  )
}

message("\n✅ All scenarios processed successfully!")
message("Final results saved to: ", output_dir)