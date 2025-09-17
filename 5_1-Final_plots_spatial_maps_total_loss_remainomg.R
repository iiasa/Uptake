# Spatial Maps with Custom Color Schemes and Column Separation
library(terra)
library(tidyverse)
library(glue)
library(sf)
library(viridis)
library(patchwork)
library(scales)
library(rnaturalearth)

# ──── CONFIGURATION ────────────────────────────────────────────────────────
scenarios <- c("Mit2p6", "Mit4p5", "Mit7p0")
decades <- c(2020, 2050, 2080, 2100)

# Directory setup
root_dir <- "H:/Uptake"
flam_dir <- "H:/UPTAKE/FLAM_Output/UKESM1-0-LL/"
output_dir <- "H:/UPTAKE/Biomass_Loss_Results_UKESM"

# European boundaries
europe <- ne_countries(scale = "medium", continent = "europe", returnclass = "sf") %>%
  st_crop(xmin = -15, xmax = 35, ymin = 35, ymax = 70)

# Conversion factors
km2_to_ha <- 100      # Convert km² to hectares (ha)
tc_to_mtc <- 1e-6     # Convert tonnes to million tonnes (MtC)

# Scenario mapping for FLAM directories
scenario_mapping <- c(
  Mit2p6 = "SSP126",
  Mit4p5 = "SSP245",
  Mit7p0 = "SSP370"
)

# ──── LOAD AND PROCESS DATA FUNCTION ──────────────────────────────────────
load_year_data <- function(year, scenario) {
  # Get corresponding SSP scenario
  ssp <- scenario_mapping[scenario]
  
  # Biomass (Total Carbon)
  biomass_file <- file.path(root_dir, scenario, glue("{scenario}_Total_C_{year}.tif"))
  if (!file.exists(biomass_file)) stop("Biomass file not found: ", biomass_file)
  biomass <- rast(biomass_file) * tc_to_mtc
  
  # Use biomass as template for extent matching
  template <- biomass
  
  # Burnt Area
  monthly_files <- list.files(
    file.path(flam_dir, paste0(ssp, "_A_burn_monthly")),
    pattern = glue("A_burn_{year}[0-9]{{2}}\\.tif$"),
    full.names = TRUE
  )
  
  if (length(monthly_files) == 0) stop("No monthly burn files found for year: ", year)
  
  # Load and resample to match biomass extent
  burn_stack <- rast(monthly_files)
  burn_stack_resampled <- resample(burn_stack, template, method = "bilinear")
  burnt_area_km2 <- sum(burn_stack_resampled, na.rm = TRUE)
  burnt_area_ha <- burnt_area_km2 * km2_to_ha
  
  # Biomass Loss
  loss_file <- file.path(output_dir, scenario, glue("Total_Loss_{year}.tif"))
  if (!file.exists(loss_file)) stop("Biomass loss file not found: ", loss_file)
  biomass_loss <- rast(loss_file) * tc_to_mtc
  biomass_loss_resampled <- resample(biomass_loss, template, method = "bilinear")
  
  # Biomass Remaining
  biomass_remaining <- biomass - biomass_loss_resampled
  
  return(list(
    biomass = biomass,
    burnt_area = burnt_area_ha,
    biomass_loss = biomass_loss_resampled,
    biomass_remaining = biomass_remaining,
    year = year
  ))
}

# ──── PLOTTING FUNCTION WITH CUSTOM COLOR SCHEMES ──────────────────────────
create_spatial_map <- function(raster, title, scale_limits, breaks, legend_name, 
                               palette = "viridis", trans = "identity") {
  # Crop to Europe
  raster_cropped <- crop(raster, europe)
  raster_masked <- mask(raster_cropped, vect(europe))
  
  # Convert to data frame
  df <- as.data.frame(raster_masked, xy = TRUE, na.rm = TRUE) %>%
    setNames(c("x", "y", "value"))
  
  # Handle transformations
  if (trans == "log10") df$value <- pmax(df$value, 0.001)
  if (trans == "sqrt") df$value <- pmax(df$value, 0)
  
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = europe, fill = NA, color = "gray30", size = 0.2) +
    scale_fill_viridis_c(
      option = palette,
      trans = trans,
      name = legend_name,
      limits = scale_limits,
      breaks = breaks,
      labels = function(x) {
        ifelse(x >= 1000, paste0(round(x/1000), "K"), format(round(x), big.mark = ","))
      },
      na.value = "transparent",
      direction = -1  # Lighter colors for lower values
    ) +
    labs(title = title) +
    coord_sf(xlim = c(-15, 35), ylim = c(35, 70)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, margin = margin(b = 3)),
      plot.margin = margin(5, 5, 15, 5)  # Extra margin at bottom for legend
    )
}

# ──── CREATE SEPARATOR LINES BETWEEN COLUMNS ──────────────────────────────
create_separator <- function() {
  ggplot() +
    geom_vline(xintercept = 0.5, color = "gray50", size = 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# ──── MAIN PROCESSING LOOP OVER SCENARIOS ─────────────────────────────────
for (scenario in scenarios) {
  # Load all data with error handling
  all_data <- list()
  for (decade in decades) {
    result <- tryCatch({
      load_year_data(decade, scenario)
    }, error = function(e) {
      message("Error loading data for ", decade, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(result)) {
      all_data[[as.character(decade)]] <- result
    }
  }
  
  # Check if we have any data
  if (length(all_data) == 0) {
    message("No data could be loaded for scenario: ", scenario)
    next
  }
  
  # Calculate GLOBAL limits for consistent color scales
  bio_vals <- c()
  loss_vals <- c()
  remaining_vals <- c()
  burnt_vals <- c()
  
  for (data in all_data) {
    bio_vals <- c(bio_vals, values(data$biomass) %>% na.omit())
    loss_vals <- c(loss_vals, values(data$biomass_loss) %>% na.omit())
    remaining_vals <- c(remaining_vals, values(data$biomass_remaining) %>% na.omit())
    burnt_vals <- c(burnt_vals, values(data$burnt_area) %>% na.omit())
  }
  
  # Set global limits for each variable
  bio_limits <- c(max(min(bio_vals, na.rm = TRUE), 1), max(bio_vals, na.rm = TRUE))
  loss_limits <- c(max(min(loss_vals, na.rm = TRUE), 0.01), max(loss_vals, na.rm = TRUE))
  remaining_limits <- c(max(min(remaining_vals, na.rm = TRUE), 1), max(remaining_vals, na.rm = TRUE))
  burnt_limits <- c(max(min(burnt_vals, na.rm = TRUE), 0.001), max(burnt_vals, na.rm = TRUE))
  
  # Calculate breaks for each variable
  bio_breaks <- pretty(bio_limits, n = 7) %>% round()
  loss_breaks <- pretty(loss_limits, n = 10) %>% round()
  remaining_breaks <- pretty(remaining_limits, n = 7) %>% round()
  burnt_breaks <- pretty(burnt_limits, n = 10) %>% round()
  
  # ──── CREATE ALL MAPS WITH CUSTOM COLOR SCHEMES ────────────────────────────
  all_plots <- list()
  
  for (i in seq_along(all_data)) {
    decade_data <- all_data[[i]]
    decade <- decade_data$year
    
    # Only add titles to top row
    titles <- if (i == 1) {
      c("Total Biomass", "Burnt Area", "Biomass Loss", "Biomass Remaining")
    } else {
      c("", "", "", "")
    }
    
    # Create maps with CUSTOM color schemes as requested
    p_biomass <- create_spatial_map(
      decade_data$biomass, titles[1], 
      scale_limits = bio_limits,
      breaks = bio_breaks,
      legend_name = "MtC",
      "viridis", "sqrt"  # DIFFERENT color for Total Biomass
    )
    
    p_burnt <- create_spatial_map(
      decade_data$burnt_area, titles[2],
      scale_limits = burnt_limits,
      breaks = burnt_breaks,
      legend_name = "ha",
      "magma", "sqrt" 
    )
    
    p_loss <- create_spatial_map(
      decade_data$biomass_loss, titles[3],
      scale_limits = loss_limits,
      breaks = loss_breaks,
      legend_name = "MtC",
      "plasma", "sqrt"
    )
    
    p_remaining <- create_spatial_map(
      decade_data$biomass_remaining, titles[4],
      scale_limits = remaining_limits,
      breaks = remaining_breaks,
      legend_name = "MtC",
      "viridis", "sqrt"  # DIFFERENT color for Biomass Remaining
    )
    
    # Create separators between columns
    sep1 <- create_separator()
    sep2 <- create_separator()
    sep3 <- create_separator()
    sep4 <- create_separator()
    
    # Add decade label on the RIGHT side
    year_label <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = decade, 
               size = 6, fontface = "bold", color = "white") +
      theme_void() +
      theme(plot.background = element_rect(fill = "black", color = NA))
    
    # Combine into single row with separators between columns
    row_plot <- p_biomass + sep1 + p_burnt + sep2 + p_loss + sep3 + p_remaining + sep4 + year_label +
      plot_layout(widths = c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 0.1))  # Added separators
    
    all_plots[[as.character(decade)]] <- row_plot
  }
  
  # ──── COMBINE ALL DECADES WITH PROPER LEGEND PLACEMENT ───────────────────
  combined_plot <- wrap_plots(all_plots, ncol = 1) +
    plot_annotation(
      title = glue("Spatial Distribution: {scenario}"),
      subtitle = "Forest Biomass, Burnt Area, Biomass Loss, and Biomass Remaining by Decade",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.spacing.x = unit(0.5, "cm"))
  
  # Save output with wider format to accommodate separators and year labels
  spatial_maps_dir <- file.path(output_dir, "spatial_maps")
  if (!dir.exists(spatial_maps_dir)) dir.create(spatial_maps_dir, recursive = TRUE)
  
  output_file <- file.path(spatial_maps_dir, glue("spatial_maps_custom_colors_separated_{scenario}.png"))
  ggsave(output_file, combined_plot, width = 24, height = 20, dpi = 300)
  
  cat("Visualization for", scenario, "saved to:", output_file, "\n")
}