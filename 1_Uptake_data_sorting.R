# ──────────────────────────────────────────────────────────────────────────────
# G4M biomass extraction 
#   • extract values
#   • converts tC ha⁻¹ → tC using landarea_ha
#   • writes yearly totals and optional rasters
# ──────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(readr)     
library(stringr)
library(purrr)
library(raster)     
library(tidyr)
library(glue)
library(sp)

# --------------------------------------------------------------------------- #
# 1. File list
# --------------------------------------------------------------------------- #
root_dir   <- "H:/Uptake"           
file_paths <- file.path(
  root_dir,
  c(
    "bioclimaDetails_t1_UPTAKE_31072025_Ref_Mit2p6_RCPref_NONE_Pco2_0.txt",
    "bioclimaDetails_t1_UPTAKE_31072025_Ref_Mit4p5_RCPref_NONE_Pco2_0.txt",
    "bioclimaDetails_t1_UPTAKE_31072025_Ref_Mit7p0_RCPref_NONE_Pco2_0.txt"
  )
)

# --------------------------------------------------------------------------- #
# 2. Helper: read one file and return a data frame of *yearly totals*
# --------------------------------------------------------------------------- #
extract_totals <- function(path) {
  
  # 2.1  read file -----------------------------------------------------------
  df <- read_delim(path, delim = "\t", col_types = cols(), trim_ws = TRUE)
  
  # 2.2  scenario label ------------------------------------------------------
  scenario <- str_extract(basename(path), "Mit[0-9]p[0-9]")
  
  # 2.3  convert per-ha → stock (t C) & aggregate ----------------------------
  df %>% 
    mutate(
      scenario = scenario,
      aff_stock_tc = biom_af_tcha * area_forest_new_ha,  # tC
      fm_stock_tc = (biom_fm_tcha * (area_forest_old_ha - (area_forest10_ha + area_forest30_ha + area_forest_primary_ha))) +
        (biom_fm10_tcha * area_forest10_ha) + 
        (biom_fm30_tcha * area_forest30_ha) + 
        (biom_fm_primary_tcha * area_forest_primary_ha)  # tC
    ) %>% 
    group_by(year, scenario) %>%  # retain scenario
    summarise(
      afforestation_tc = sum(aff_stock_tc, na.rm = TRUE),
      forest_mgmt_tc = sum(fm_stock_tc, na.rm = TRUE),
      total_tc = afforestation_tc + forest_mgmt_tc,
      .groups = "drop"
    )
}

# --------------------------------------------------------------------------- #
# 3. Run for all scenarios and export ---------------------------------------- #
totals_tc <- map_dfr(file_paths, extract_totals)

totals_gtc <- totals_tc %>% 
  mutate(across(ends_with("_tc"), ~ .x / 1e9,  # Convert to Gt C
                .names = "{str_remove(.col, '_tc')}_GtC"))

write_csv(totals_tc, file.path(root_dir, "G4M_biomass_totals_tc.csv"))
write_csv(totals_gtc, file.path(root_dir, "G4M_biomass_totals_GtC.csv"))

message("✓ Yearly totals written to CSV")

# --------------------------------------------------------------------------- #
# 4. Write rasters of *per-ha* biomass 
# --------------------------------------------------------------------------- #
write_rasters <- TRUE 
if (write_rasters) {
  
  template_raster <- raster(
    extent(-25, 45, 35, 70),  
    res = 0.1,
    crs = "+proj=longlat +datum=WGS84 +no_defs"
  )
  
  process_for_raster <- function(path) {
    
    df <- read_delim(path, delim = "\t", col_types = cols(), trim_ws = TRUE)
    scenario <- str_extract(basename(path), "Mit[0-9]p[0-9]")
    out_dir <- file.path(root_dir, scenario)
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    years <- unique(df$year)
    
    for (yr in years) {
      yearly <- df %>% 
        filter(year == yr) %>% 
        mutate(
          Biom_af = biom_af_tcha * area_forest_new_ha,  # tC 
          Biom_fm_total = (biom_fm_tcha * (area_forest_old_ha - (area_forest10_ha + area_forest30_ha + area_forest_primary_ha))) +
            (biom_fm10_tcha * area_forest10_ha) + 
            (biom_fm30_tcha * area_forest30_ha) + 
            (biom_fm_primary_tcha * area_forest_primary_ha),
          Deadwood_old_total = deadwood_old_tCha * area_forest_old_ha,
          Deadwood_new_total = deadwood_new_tCha * area_forest_new_ha,
          Litter_old_total = litter_old_tCha * area_forest_old_ha,
          Litter_new_total = litter_new_tCha * area_forest_new_ha
        )
      
      coords <- yearly %>% dplyr::select(x, y)
      sp <- SpatialPointsDataFrame(coords, yearly, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      # Create and write rasters
      raster_fields <- c("Biom_af", "Biom_fm_total", 
                         "Deadwood_old_total", "Deadwood_new_total",
                         "Litter_old_total", "Litter_new_total")
      
      for (field in raster_fields) {
        r <- rasterize(sp, template_raster, field = field, fun = sum)
        writeRaster(
          r,
          filename = file.path(out_dir, glue("{scenario}_{field}_{yr}.tif")),
          format = "GTiff", 
          overwrite = TRUE
        )
      }
    }
    message("✓ Rasters written for ", scenario)
  }
  
  walk(file_paths, process_for_raster)
}