# ──────────────────────────────────────────────────────────────────────────────
# Comprehensive Biomass Analysis: Stocks, Loss, Burnt Area & Remaining Biomass
#   • Calculates total biomass (existing + afforested) in MtC
#   • Computes biomass loss from wildfires using main/shadow FLAM inputs
#   • Determines remaining biomass after losses
#   • Generates time series plots and spatial maps for key years
#   • Exports all results as CSV files
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
  library(glue)
  library(rnaturalearth)
  library(sf)
  library(viridis)
  library(patchwork)
  library(scales)
})

# ──── CONFIGURATION ────────────────────────────────────────────────────────
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
flam_scenarios <- c("SSP126", "SSP245", "SSP370")
years <- 1990:2100
map_years <- c(2020, 2050, 2080, 2100)

# Directory setup
g4m_dir <- "H:/Uptake"
prop_dir <- "H:/UPTAKE/Area_Proportions"
output_dir <- "H:/UPTAKE/Biomass_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# FLAM directories (main and shadow/uncertainty)
flam_dirs <- c(
  main = "H:/UPTAKE/FLAM_Output/UKESM1-0-LL",
  shadow = "H:/UPTAKE/FLAM_Output"  # Shadow/uncertainty runs
)

# European boundaries
europe <- ne_countries(scale = "medium", continent = "europe", returnclass = "sf") %>%
  st_crop(xmin = -15, xmax = 35, ymin = 35, ymax = 70)

# Create scenario mapping
scenario_map <- tibble(
  g4m_scenario = scenarios,
  flam_scenario = flam_scenarios
)

# Template raster
template <- rast(ext(-25, 45, 35, 70), res = 0.5, crs = "EPSG:4326")

# ──── BIOMASS CALCULATION FUNCTIONS ────────────────────────────────────────

# Calculate biomass stocks in MtC (1e6 tonnes)
calculate_biomass_stocks <- function(scenario) {
  message("Processing biomass stocks for: ", scenario)
  
  # Create output directories upfront
  scenario_dir <- file.path(output_dir, scenario)
  map_dir <- file.path(scenario_dir, "Biomass_Maps")
  dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)
  
  biomass_data <- tibble()
  
  for (yr in years) {
    # Load biomass rasters
    aff_file <- file.path(g4m_dir, scenario, glue("{scenario}_Biom_af_{yr}.tif"))
    mgd_file <- file.path(g4m_dir, scenario, glue("{scenario}_Biom_fm_total_{yr}.tif"))
    
    if (!file.exists(aff_file) || !file.exists(mgd_file)) {
      message("Skipping year ", yr, " - files not found")
      next
    }
    
    tryCatch({
      r_aff <- rast(aff_file)
      r_mgd <- rast(mgd_file)
      
      # Convert to MtC (1e6 tonnes)
      aff_sum <- global(r_aff, sum, na.rm = TRUE)[1,1] / 1e6
      mgd_sum <- global(r_mgd, sum, na.rm = TRUE)[1,1] / 1e6
      total_sum <- aff_sum + mgd_sum
      
      biomass_data <- rbind(biomass_data, tibble(
        year = yr,
        scenario = scenario,
        total_biomass = total_sum,
        aff_biomass = aff_sum,
        mgd_biomass = mgd_sum
      ))
      
      # Save spatial maps for target years
      if (yr %in% map_years) {
        total_biomass <- r_aff + r_mgd
        
        # Save raster
        writeRaster(
          total_biomass,
          file.path(map_dir, glue("Total_Biomass_{yr}.tif")),
          overwrite = TRUE
        )
        
        # Create map plot
        biomass_df <- as.data.frame(total_biomass, xy = TRUE) %>%
          setNames(c("x", "y", "biomass"))
        
        p <- ggplot() +
          geom_tile(data = biomass_df, aes(x, y, fill = biomass/1e6)) + # Convert to MtC
          geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
          scale_fill_viridis_c(
            option = "viridis",
            trans = "log10",
            name = "Biomass (MtC)",
            labels = scales::comma,
            na.value = NA
          ) +
          coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
          labs(title = glue("Total Biomass: {scenario} {yr}"),
               x = "", y = "") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(
          file.path(map_dir, glue("Total_Biomass_Map_{yr}.png")),
          p,
          width = 10,
          height = 8,
          dpi = 300
        )
      }
    }, error = function(e) {
      message("Error processing year ", yr, ": ", e$message)
    })
  }
  
  # Save time series data
  write_csv(
    biomass_data,
    file.path(scenario_dir, "Biomass_Stocks.csv")
  )
  
  # Create time series plot
  if (nrow(biomass_data) > 0) {
    p_time <- biomass_data %>%
      pivot_longer(
        cols = -c(year, scenario),
        names_to = "type",
        values_to = "biomass"
      ) %>%
      mutate(type = factor(
        type,
        levels = c("total_biomass", "mgd_biomass", "aff_biomass"),
        labels = c("Total", "Managed", "Afforested")
      )) %>%
      ggplot(aes(year, biomass, color = type)) +
      geom_line(size = 1.2) +
      geom_point(size = 2) +
      scale_color_manual(
        values = c("Total" = "#1f77b4", "Managed" = "#ff7f0e", "Afforested" = "#2ca02c")
      ) +
      labs(
        title = glue("Biomass Stocks: {scenario}"),
        x = "Year",
        y = "Biomass (MtC)",
        color = "Forest Type"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90"),
        plot.title = element_text(hjust = 0.5)
      )
    
    ggsave(
      file.path(scenario_dir, "Biomass_Stocks_TimeSeries.png"),
      p_time,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  return(biomass_data)
}

# Update the biomass loss calculation function
calculate_biomass_loss <- function(g4m_scenario, flam_scenario, flam_source) {
  scenario_out_dir <- file.path(output_dir, g4m_scenario, "Biomass_Loss", flam_source)
  dir.create(scenario_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  annual_results <- tibble()
  spatial_data <- list()
  
  for (year in 2015:2100) {
    tryCatch({
      # 1. FLAM BURNT AREA PROCESSING ───────────────────────────────────────
      flam_path <- file.path(flam_dirs[flam_source], glue("{flam_scenario}_A_burn_monthly"))
      
      if (!dir.exists(flam_path)) {
        message("FLAM directory not found: ", flam_path)
        next
      }
      
      flam_files <- list.files(
        flam_path,
        pattern = glue("{year}.*\\.tif$"),
        full.names = TRUE
      )
      
      if (length(flam_files) == 0) {
        message("No FLAM files found for ", year, " in ", flam_path)
        next
      }
      
      # Calculate annual burnt area
      burn_stack <- rast(flam_files)
      annual_burn <- sum(burn_stack)
      
      # Resample FLAM data to template grid
      annual_burn_resampled <- resample(annual_burn, template, method = "bilinear")
      
      # 2. FOREST PROPORTIONS ───────────────────────────────────────────────
      aff_prop <- rast(file.path(prop_dir, g4m_scenario,
                                 glue("{g4m_scenario}_afforested_proportion_stack.tif")))
      mgd_prop <- rast(file.path(prop_dir, g4m_scenario,
                                 glue("{g4m_scenario}_managed_proportion_stack.tif")))
      
      aff_prop_yr <- aff_prop[[glue("y{year}")]]
      mgd_prop_yr <- mgd_prop[[glue("y{year}")]]
      
      # 3. BURNT AREA BY FOREST TYPE ────────────────────────────────────────
      aff_burn <- annual_burn_resampled * aff_prop_yr
      mgd_burn <- annual_burn_resampled * mgd_prop_yr
      
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
        annual_burn_resampled,
        file.path(scenario_out_dir, glue("Total_Burn_{year}.tif")),
        overwrite = TRUE
      )
      
      # 4. BIOMASS DENSITY PROCESSING ───────────────────────────────────────
      aff_dens <- rast(file.path(g4m_dir, g4m_scenario, glue("{g4m_scenario}_Biom_af_{year}.tif")))
      mgd_dens <- rast(file.path(g4m_dir, g4m_scenario, glue("{g4m_scenario}_Biom_fm_total_{year}.tif")))
      
      # Resample biomass to template grid
      aff_dens_resampled <- resample(aff_dens, template, method = "bilinear")
      mgd_dens_resampled <- resample(mgd_dens, template, method = "bilinear")
      
      # Convert to tC/ha using template cell size
      cell_area_ha <- cellSize(template, unit = "ha")
      aff_dens_tcha <- aff_dens_resampled / cell_area_ha
      mgd_dens_tcha <- mgd_dens_resampled / cell_area_ha
      
      # 5. BIOMASS LOSS CALCULATION ────────────────────────────────────────
      aff_loss <- aff_burn * aff_dens_tcha * 100  # Convert to tonnes
      mgd_loss <- mgd_burn * mgd_dens_tcha * 100
      total_loss <- aff_loss + mgd_loss
      
      # 6. REMAINING BIOMASS CALCULATION ──────────────────────────────────
      total_biomass <- aff_dens_resampled + mgd_dens_resampled  # Total biomass in tonnes
      remaining_aff <- aff_dens_resampled - aff_loss
      remaining_mgd <- mgd_dens_resampled - mgd_loss
      remaining_total <- total_biomass - total_loss
      
      # Convert to MtC for summaries
      aff_loss_mtc <- sum(values(aff_loss), na.rm = TRUE) / 1e6
      mgd_loss_mtc <- sum(values(mgd_loss), na.rm = TRUE) / 1e6
      total_loss_mtc <- sum(values(total_loss), na.rm = TRUE) / 1e6
      aff_remain_mtc <- sum(values(remaining_aff), na.rm = TRUE) / 1e6
      mgd_remain_mtc <- sum(values(remaining_mgd), na.rm = TRUE) / 1e6
      total_remain_mtc <- sum(values(remaining_total), na.rm = TRUE) / 1e6
      
      # Store results
      annual_results <- rbind(annual_results, tibble(
        year = year,
        scenario = g4m_scenario,
        source = flam_source,
        total_burn = sum(values(annual_burn_resampled), na.rm = TRUE),
        aff_burn = sum(values(aff_burn), na.rm = TRUE),
        mgd_burn = sum(values(mgd_burn), na.rm = TRUE),
        aff_loss = aff_loss_mtc,
        mgd_loss = mgd_loss_mtc,
        total_loss = total_loss_mtc,
        aff_remain = aff_remain_mtc,
        mgd_remain = mgd_remain_mtc,
        total_remain = total_remain_mtc
      ))
      
      # Save spatial data for key years
      if (year %in% map_years) {
        spatial_data[[as.character(year)]] <- list(
          total_loss = total_loss,
          total_remain = remaining_total
        )
        
        # Save remaining biomass raster
        writeRaster(
          remaining_total,
          file.path(scenario_out_dir, glue("Remaining_Biomass_{year}.tif")),
          overwrite = TRUE
        )
      }
    }, error = function(e) {
      message("Error processing year ", year, ": ", e$message)
    })
  }
  
  # Save annual results
  write_csv(annual_results, file.path(scenario_out_dir, "annual_results.csv"))
  return(annual_results)
}
#Function to create spatial comparison maps
create_spatial_comparison <- function(scenario, year, source_type) {
  # Load biomass data
  biomass_file <- file.path(output_dir, scenario, "Biomass_Maps", 
                            glue("Total_Biomass_{year}.tif"))
  if (!file.exists(biomass_file)) return(NULL)
  
  biomass <- rast(biomass_file)
  
  # Load remaining biomass
  remain_file <- file.path(output_dir, scenario, "Biomass_Loss", source_type,
                           glue("Remaining_Biomass_{year}.tif"))
  if (!file.exists(remain_file)) return(NULL)
  
  remaining <- rast(remain_file)
  
  # Resample both to template grid
  biomass_resampled <- resample(biomass, template, method = "bilinear")
  remaining_resampled <- resample(remaining, template, method = "bilinear")
  
  # Create ocean mask
  europe_vect <- vect(europe)
  ocean_mask <- rasterize(europe_vect, template, background = NA)
  
  # Apply masks
  biomass_masked <- mask(biomass_resampled, ocean_mask, inverse = FALSE)
  remain_masked <- mask(remaining_resampled, ocean_mask, inverse = FALSE)
  
  # Convert to MtC
  biomass_mtc <- biomass_masked / 1e6
  remain_mtc <- remain_masked / 1e6
  
  # Create data frames
  biomass_df <- as.data.frame(biomass_mtc, xy = TRUE, na.rm = FALSE) %>%
    setNames(c("x", "y", "biomass"))
  remain_df <- as.data.frame(remain_mtc, xy = TRUE, na.rm = FALSE) %>%
    setNames(c("x", "y", "remaining"))
  
  # Create plots
  p_biomass <- ggplot() +
    geom_tile(data = biomass_df, aes(x, y, fill = biomass)) +
    geom_sf(data = europe, fill = NA, color = "gray30", size = 0.2) +
    scale_fill_viridis_c(
      option = "viridis",
      name = "Biomass (MtC)",
      na.value = "transparent"
    ) +
    coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
    labs(title = glue("Total Biomass: {scenario} {year}")) +
    theme_void()
  
  p_remain <- ggplot() +
    geom_tile(data = remain_df, aes(x, y, fill = remaining)) +
    geom_sf(data = europe, fill = NA, color = "gray30", size = 0.2) +
    scale_fill_viridis_c(
      option = "viridis",
      name = "Remaining Biomass (MtC)",
      na.value = "transparent"
    ) +
    coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
    labs(title = glue("Remaining Biomass: {scenario} {year} ({source_type})")) +
    theme_void()
  
  # Combine plots
  p_combined <- p_biomass + p_remain +
    plot_annotation(title = glue("Biomass Comparison: {scenario} {year}"),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
  
  # Save output
  out_dir <- file.path(output_dir, scenario, "Spatial_Comparisons")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    file.path(out_dir, glue("Biomass_Comparison_{year}_{source_type}.png")),
    p_combined,
    width = 16,
    height = 8,
    dpi = 300
  )
}

# ──── MAIN PROCESSING ─────────────────────────────────────────────────────

# 1. Calculate biomass stocks for all scenarios
all_biomass <- map_dfr(scenarios, calculate_biomass_stocks)

# 2. Calculate biomass loss and remaining biomass for all scenarios/sources
all_loss <- map_dfr(1:nrow(scenario_map), function(i) {
  g4m_scen <- scenario_map$g4m_scenario[i]
  flam_scen <- scenario_map$flam_scenario[i]
  
  map_dfr(names(flam_dirs), function(src) {
    calculate_biomass_loss(g4m_scen, flam_scen, src)
  })
})

# 3. Generate comparison plots ─────────────────────────────────────────────

# Biomass stocks comparison
p_stocks <- ggplot(all_biomass, aes(year, total_biomass, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c("Mit2p6" = "#1b9e77", "Mit4p5" = "#d95f02", "Mit7p0" = "#7570b3")
  ) +
  labs(
    title = "Total Biomass Comparison",
    x = "Year",
    y = "Biomass (MtC)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)

# Biomass loss comparison (main source only)
p_loss <- all_loss %>%
  filter(source == "main") %>%
  ggplot(aes(year, total_loss, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c("Mit2p6" = "#1b9e77", "Mit4p5" = "#d95f02", "Mit7p0" = "#7570b3")
  ) +
  labs(
    title = "Biomass Loss from Wildfires",
    x = "Year",
    y = "Biomass Loss (MtC)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)

# Remaining biomass comparison
p_remain <- all_loss %>%
  filter(source == "main") %>%
  ggplot(aes(year, total_remain, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c("Mit2p6" = "#1b9e77", "Mit4p5" = "#d95f02", "Mit7p0" = "#7570b3")
  ) +
  labs(
    title = "Remaining Biomass After Wildfires",
    x = "Year",
    y = "Biomass (MtC)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)

# Burnt area comparison
p_burn <- all_loss %>%
  filter(source == "main") %>%
  ggplot(aes(year, total_burn, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c("Mit2p6" = "#1b9e77", "Mit4p5" = "#d95f02", "Mit7p0" = "#7570b3")
  ) +
  labs(
    title = "Total Burnt Area",
    x = "Year",
    y = "Burnt Area (km²)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)

# Save comparison plots
walk2(
  list(p_stocks, p_loss, p_remain, p_burn),
  c("Stocks_Comparison", "Loss_Comparison", "Remain_Comparison", "Burn_Comparison"),
  ~ ggsave(
    file.path(output_dir, glue("{.y}.png")),
    .x,
    width = 10,
    height = 6,
    dpi = 300
  )
)

# 4. Generate spatial comparison maps ─────────────────────────────────────
for (scen in scenarios) {
  for (yr in map_years) {
    for (src in names(flam_dirs)) {
      create_spatial_comparison(scen, yr, src)
    }
  }
}

# 5. Export all data to CSV ───────────────────────────────────────────────

# Biomass stocks
write_csv(all_biomass, file.path(output_dir, "All_Biomass_Stocks.csv"))

# Biomass loss and remaining biomass
write_csv(all_loss, file.path(output_dir, "All_Biomass_Loss.csv"))

message("\n✅ All processing completed successfully!")
message("Results saved to: ", output_dir)




























# ──────────────────────────────────────────────────────────────────────────────
# Final Enhanced Cumulative Biomass Analysis with Data Export
#   • Added CSV export for all calculated data
#   • Professional styling with boxed plots
#   • Comprehensive data output for further analysis
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(glue)
library(patchwork)
library(scales)
library(ggthemes)
library(ggtext)

# ──── CONFIGURATION ────────────────────────────────────────────────────────
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
years <- 2015:2100
output_dir <- "H:/UPTAKE/Biomass_Results"
adjustment_factor <- 25  # 0.5° to 0.1° resolution (0.5/0.1)^2 = 25

# Enhanced color palette
palette <- c(
  "Afforested" = "#1b9e77",  # Green
  "Managed" = "#d95f02"       # Orange
)

# Font settings
font_title <- "sans"
font_base <- "sans"

# ──── LOAD AND PROCESS DATA ────────────────────────────────────────────────

# Function to load and process scenario data
load_scenario_data <- function(scenario) {
  # Load biomass stocks
  stocks_file <- file.path(output_dir, scenario, "Biomass_Stocks.csv")
  stocks <- read_csv(stocks_file) %>%
    dplyr::select(year, aff_biomass, mgd_biomass) %>%
    filter(year >= 2015) %>%
    mutate(
      aff_biomass = aff_biomass,  # Already in MtC
      mgd_biomass = mgd_biomass
    )
  
  # Load loss data for both sources
  main_loss <- read_csv(file.path(output_dir, scenario, "Biomass_Loss/main/annual_results.csv"))
  shadow_loss <- read_csv(file.path(output_dir, scenario, "Biomass_Loss/shadow/annual_results.csv"))
  
  # Combine and calculate cumulative values
  cumulative_data <- bind_rows(
    main_loss %>% mutate(source = "main"),
    shadow_loss %>% mutate(source = "shadow")
  ) %>%
    arrange(year) %>%
    # Join with stocks data
    left_join(stocks, by = "year") %>%
    group_by(source) %>%
    mutate(
      # Apply resolution adjustment
      across(c(ends_with("burn"), ends_with("loss"), ends_with("remain")), 
             ~ . * adjustment_factor),
      
      # Apply additional scaling corrections
      across(c(ends_with("burn")), ~ . / 100),         # Convert to Mha (1 Mha = 100 km²)
      across(c(ends_with("loss")), ~ . * 10000),        # Correct 10000x underestimation
      across(c(ends_with("loss"), ends_with("remain")), ~ . / 1e6),  # Convert to MtC
      
      # Calculate cumulative values
      aff_cum_burn = cumsum(aff_burn),
      mgd_cum_burn = cumsum(mgd_burn),
      
      aff_cum_loss = cumsum(aff_loss),
      mgd_cum_loss = cumsum(mgd_loss),
      
      # Remaining biomass (stocks - cumulative loss)
      aff_remain = aff_biomass - aff_cum_loss,
      mgd_remain = mgd_biomass - mgd_cum_loss
    ) %>%
    ungroup()
  
  return(cumulative_data)
}

# Load all data
all_data <- map_dfr(scenarios, ~ {
  data <- load_scenario_data(.x)
  data$scenario <- .x
  return(data)
})

# ──── EXPORT COMBINED DATA TO CSV ──────────────────────────────────────────

# Create output directory
data_dir <- file.path(output_dir, "Combined_Data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Function to prepare and export component data
export_component_data <- function(component, file_name) {
  # Prepare data
  component_data <- all_data %>%
    dplyr::select(year, scenario, source,
                  aff_value = !!sym(paste0("aff_", component)),
                  mgd_value = !!sym(paste0("mgd_", component))) %>%
    pivot_longer(
      cols = c(aff_value, mgd_value),
      names_to = "forest_type",
      values_to = "value"
    ) %>%
    mutate(
      forest_type = case_when(
        forest_type == "aff_value" ~ "Afforested",
        forest_type == "mgd_value" ~ "Managed"
      ),
      component = component
    )
  
  # Save to CSV
  write_csv(component_data, file.path(data_dir, file_name))
  
  return(component_data)
}

# Export all components
biomass_data <- export_component_data("biomass", "Biomass_Stocks.csv")
burn_data <- export_component_data("cum_burn", "Cumulative_Burnt_Area.csv")
loss_data <- export_component_data("cum_loss", "Cumulative_Biomass_Loss.csv")
remain_data <- export_component_data("remain", "Remaining_Biomass.csv")

# Combine all data for comprehensive export
all_export_data <- bind_rows(
  biomass_data,
  burn_data,
  loss_data,
  remain_data
) %>%
  mutate(
    component = case_when(
      component == "biomass" ~ "Biomass Stocks",
      component == "cum_burn" ~ "Cumulative Burnt Area",
      component == "cum_loss" ~ "Cumulative Biomass Loss",
      component == "remain" ~ "Remaining Biomass",
      TRUE ~ component
    )
  )

# Save comprehensive dataset
write_csv(all_export_data, file.path(data_dir, "All_Cumulative_Data.csv"))

# ──── CREATE PROFESSIONAL DUAL-PANEL PLOTS ────────────────────────────────

# Function to create professional dual-panel component plot
create_professional_plot <- function(component, y_label, title) {
  # Prepare data
  plot_data <- all_data %>%
    dplyr::select(year, scenario, 
                  aff_value = !!sym(paste0("aff_", component)),
                  mgd_value = !!sym(paste0("mgd_", component))) %>%
    pivot_longer(
      cols = c(aff_value, mgd_value),
      names_to = "forest_type",
      values_to = "value"
    ) %>%
    mutate(
      forest_type = case_when(
        forest_type == "aff_value" ~ "Afforested",
        forest_type == "mgd_value" ~ "Managed"
      ),
      # Create factor for ordering
      forest_type = factor(forest_type, levels = c("Managed", "Afforested"))
    )
  
  # Create uncertainty data
  uncertainty <- plot_data %>%
    group_by(year, scenario, forest_type) %>%
    summarise(
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Remove infinite values
    mutate(
      min_val = ifelse(is.infinite(min_val), NA, min_val),
      max_val = ifelse(is.infinite(max_val), NA, max_val)
    )
  
  # Create plot
  p <- ggplot(uncertainty, aes(x = year, color = forest_type, fill = forest_type)) +
    geom_ribbon(aes(ymin = min_val, ymax = max_val), alpha = 0.2) +
    geom_line(aes(y = mean_val), size = 1.2) +
    facet_grid(forest_type ~ scenario, scales = "free_y") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(
      title = title,
      x = "Year",
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "black", size = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        family = font_title, 
        hjust = 0.5, 
        face = "bold",
        size = 16,
        margin = margin(b = 10)
      ),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(family = font_base, face = "bold"),
      axis.text = element_text(color = "black"),
      plot.margin = margin(15, 15, 15, 15),
      strip.background = element_rect(fill = "gray90", color = "black", size = 0.5)
    ) +
    scale_x_continuous(breaks = seq(2020, 2100, by = 20))
  
  return(p)
}

# Create individual component plots with professional styling
p_stocks <- create_professional_plot("biomass", "Biomass (MtC)", "Biomass Stocks")
p_burn <- create_professional_plot("cum_burn", "Area (Mha)", "Cumulative Burnt Area")
p_loss <- create_professional_plot("cum_loss", "Loss (MtC)", "Cumulative Biomass Loss")
p_remain <- create_professional_plot("remain", "Biomass (MtC)", "Remaining Biomass")

# ──── CREATE ENHANCED COMBINED PLOT ────────────────────────────────────────

# Combine all plots in a 2x2 grid
p_combined <- (p_stocks + p_burn) / (p_loss + p_remain) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

# Add global title and annotations
p_combined <- p_combined +
  plot_annotation(
    title = "Comprehensive Biomass Analysis Across Climate Scenarios",
    subtitle = "Showing cumulative metrics with uncertainty ranges (2015-2100)",
    caption = "Data Sources: FLAM-G4M Integration",
    theme = theme(
      plot.title = element_text(
        family = font_title,
        hjust = 0.5,
        size = 24,
        face = "bold",
        margin = margin(t = 10, b = 5)
      ),
      plot.subtitle = element_text(
        family = font_base,
        hjust = 0.5,
        size = 18,
        color = "gray30",
        margin = margin(b = 15)
      ),
      plot.caption = element_text(
        family = font_base,
        hjust = 1,
        size = 12,
        color = "gray40",
        margin = margin(t = 10)
      ),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  )

# Save output
out_dir <- file.path(output_dir, "Enhanced_Analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  file.path(out_dir, "Enhanced_Cumulative_Analysis.png"),
  p_combined,
  width = 16,
  height = 14,
  dpi = 600
)

# ──── CREATE SCENARIO-SPECIFIC PLOTS ────────────────────────────────────────

# Function to create professional scenario-specific plot
create_scenario_plot <- function(scen) {
  # Filter data for scenario
  scenario_data <- all_data %>% filter(scenario == scen)
  
  # Create plot data
  plot_data <- scenario_data %>%
    dplyr::select(year, 
                  aff_biomass, mgd_biomass,
                  aff_cum_burn, mgd_cum_burn,
                  aff_cum_loss, mgd_cum_loss,
                  aff_remain, mgd_remain) %>%
    pivot_longer(
      cols = -year,
      names_to = c("forest_type", "metric"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    mutate(
      forest_type = case_when(
        forest_type == "aff" ~ "Afforested",
        forest_type == "mgd" ~ "Managed"
      ),
      metric = case_when(
        metric == "biomass" ~ "Biomass Stocks",
        metric == "cum_burn" ~ "Cumulative Burnt Area",
        metric == "cum_loss" ~ "Cumulative Biomass Loss",
        metric == "remain" ~ "Remaining Biomass",
        TRUE ~ metric
      ),
      # Create factor for ordering
      forest_type = factor(forest_type, levels = c("Managed", "Afforested"))
    )
  
  # Create uncertainty data
  uncertainty <- plot_data %>%
    group_by(year, forest_type, metric) %>%
    summarise(
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Remove infinite values
    mutate(
      min_val = ifelse(is.infinite(min_val), NA, min_val),
      max_val = ifelse(is.infinite(max_val), NA, max_val)
    )
  
  # Create plot
  p <- ggplot(uncertainty, aes(x = year, color = forest_type, fill = forest_type)) +
    geom_ribbon(aes(ymin = min_val, ymax = max_val), alpha = 0.2) +
    geom_line(aes(y = mean_val), size = 1.2) +
    facet_grid(forest_type ~ metric, scales = "free_y") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(
      title = glue("Cumulative Biomass Analysis: {scen} Scenario"),
      subtitle = "With uncertainty from FLAM sources",
      x = "Year",
      y = "Value"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "black", size = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        family = font_title,
        hjust = 0.5,
        face = "bold",
        size = 20,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        family = font_base,
        hjust = 0.5,
        size = 16,
        color = "gray30",
        margin = margin(b = 15)
      ),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "gray90", color = "black", size = 0.5),
      axis.title = element_text(family = font_base, face = "bold"),
      axis.text = element_text(color = "black"),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_x_continuous(breaks = seq(2020, 2100, by = 20))
  
  ggsave(
    file.path(out_dir, glue("Enhanced_Analysis_{scen}.png")),
    p,
    width = 16,
    height = 10,
    dpi = 600
  )
}

# Create scenario-specific plots
for (scen in scenarios) {
  create_scenario_plot(scen)
}

message("✅ All enhanced plots created successfully!")
message("✅ All data exported to CSV files")
message("Results saved to: ", out_dir)