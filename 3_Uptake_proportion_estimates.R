# ──────────────────────────────────────────────────────────────────────────────
# Area Proportions with Forest Mask + Renormalization
#   • Loads G4M table once per scenario
#   • Builds forest mask per year (total_forest_ha > 0)
#   • Resamples with nearest neighbour (no bilinear dilution)
#   • Renormalizes aff + mgd = 1 inside forest
#   • Saves stacks and runs forest-only sanity checks
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(terra)
  library(glue)
  library(ggplot2)
})

# ------------------------------- Config ------------------------------------- #
root_dir        <- "H:/Uptake"
root_dir_output <- "H:/Uptake/Area_Proportions"
scenarios       <- c("Mit2p6","Mit4p5","Mit7p0")          # add more if needed
years           <- 1990:2100
tolerance       <- 1e-3                  # sanity tolerance for aff+mgd
template        <- rast(ext(-25, 45, 35, 70), res = 0.5, crs = "EPSG:4326")

if (!dir.exists(root_dir_output)) dir.create(root_dir_output, recursive = TRUE)

# -------------------------- Helper Functions -------------------------------- #
create_xyz_rast <- function(data, field) {
  # data has columns x, y, <field>
  rast(data[, c("x", "y", field)], type = "xyz", crs = "EPSG:4326")
}

resample_near <- function(r, tmpl) {
  # nearest neighbour to avoid mixing zeros and forest
  resample(r, tmpl, method = "near")
}

save_stack <- function(stack_list, years_vec, out_dir, scenario, name) {
  if (length(stack_list) == 0) return(invisible(NULL))
  s <- rast(stack_list)
  # use actual year keys to avoid misalignment if some years skipped
  names(s) <- paste0("y", years_vec)
  writeRaster(
    s,
    file.path(out_dir, glue("{scenario}_{name}_proportion_stack.tif")),
    overwrite = TRUE,
    filetype  = "GTiff",
    datatype  = "FLT4S",
    NAflag    = -9999
  )
}

# ------------------------------ Main Loop ----------------------------------- #
for (scenario in scenarios) {
  message("\n────────────────────────────────────────────────")
  message("Processing scenario: ", scenario)

  # Input table (once)
  input_file <- file.path(
    root_dir,
    glue("bioclimaDetails_t1_UPTAKE_31072025_Ref_{scenario}_RCPref_NONE_Pco2_0.txt")
  )
  if (!file.exists(input_file)) {
    warning("Missing input for scenario: ", scenario)
    next
  }

  df_all <- suppressMessages(read_delim(input_file, delim = "\t", trim_ws = TRUE))
  needed_cols <- c("x","y","year","area_forest_old_ha","area_forest_new_ha")
  miss <- setdiff(needed_cols, names(df_all))
  if (length(miss)) stop("Missing columns in input: ", paste(miss, collapse = ", "))

  # Containers
  aff_stack <- list()
  mgd_stack <- list()
  yrs_done  <- c()

  # Per-year processing
  for (yr in years) {
    df_year <- df_all %>% filter(year == !!yr)
    if (nrow(df_year) == 0) {
      message("Skipping year (no rows): ", yr)
      next
    }

    df_props <- df_year %>%
      mutate(
        total_forest_ha = area_forest_old_ha + area_forest_new_ha,
        aff_prop_raw    = if_else(total_forest_ha > 0,
                                  area_forest_new_ha / total_forest_ha, 0),
        mgd_prop_raw    = if_else(total_forest_ha > 0,
                                  area_forest_old_ha / total_forest_ha, 0)
      ) %>%
      mutate(
        aff_prop_raw = pmax(0, pmin(1, aff_prop_raw)),
        mgd_prop_raw = pmax(0, pmin(1, mgd_prop_raw)),
        forest_flag  = as.integer(total_forest_ha > 0)
      )

    # Native rasters
    r_aff_native   <- create_xyz_rast(df_props, "aff_prop_raw")
    r_mgd_native   <- create_xyz_rast(df_props, "mgd_prop_raw")
    r_mask_native  <- create_xyz_rast(df_props, "forest_flag")

    # Resample with nearest neighbour (no mixing)
    r_aff <- resample_near(r_aff_native,  template)
    r_mgd <- resample_near(r_mgd_native,  template)
    r_msk <- resample_near(r_mask_native, template)

    # Renormalize inside forest so that aff + mgd = 1
    # Keep zeros outside forest (msk == 0); set NA for strict "non-forest"
    s <- r_aff + r_mgd
    inside <- r_msk == 1 & s > 0
    # avoid division by zero; only renormalize where s>0 & forest
    r_aff_norm <- ifel(inside, r_aff / s, 0)
    r_mgd_norm <- ifel(inside, r_mgd / s, 0)

    # Optional—strict NA outside forest for analysis (kept zeros in saved stacks)
    # r_aff_norm <- mask(r_aff_norm, r_msk, maskvalues = 0)
    # r_mgd_norm <- mask(r_mgd_norm, r_msk, maskvalues = 0)

    aff_stack[[as.character(yr)]] <- r_aff_norm
    mgd_stack[[as.character(yr)]] <- r_mgd_norm
    yrs_done  <- c(yrs_done, yr)

    if (yr %% 10 == 0) message(" …year ", yr, " done")
  }

  if (!length(yrs_done)) next

  # Save stacks
  out_dir <- file.path(root_dir_output, scenario)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_stack(aff_stack, yrs_done, out_dir, scenario, "afforested")
  save_stack(mgd_stack, yrs_done, out_dir, scenario, "managed")

  # ---------------------- Forest-only Sanity Checks ------------------------- #
  message("Running forest-only sanity checks…")

  aff_stack_r <- rast(file.path(out_dir, glue("{scenario}_afforested_proportion_stack.tif")))
  mgd_stack_r <- rast(file.path(out_dir, glue("{scenario}_managed_proportion_stack.tif")))
  stopifnot(nlyr(aff_stack_r) == nlyr(mgd_stack_r))

  # Build a yearly forest mask from the stacks (sum>0 means forest)
  sum_stack <- aff_stack_r + mgd_stack_r
  # mean within-forest; ignore non-forest by masking NAs for sum==0
  sum_stack_masked <- classify(sum_stack, rcl = matrix(c(-Inf, 0, NA), ncol = 3, byrow = TRUE), right = TRUE)

  mean_sum <- global(sum_stack_masked, fun = "mean", na.rm = TRUE)[,1]
  yrs_lab  <- as.numeric(sub("^y", "", names(sum_stack)))

  # Fraction invalid inside forest per year
  # invalid = |aff+mgd-1| > tolerance
  invalid_stack <- abs(sum_stack - 1) > tolerance
  # mask to forest
  invalid_inside <- mask(invalid_stack, sum_stack, maskvalues = 0)
  # share of invalid among forest cells:
  # (count TRUE) / (count non-NA)
  count_true   <- global(invalid_inside, fun = "sum", na.rm = TRUE)[,1]
  count_forest <- global(!is.na(sum_stack_masked), fun = "sum", na.rm = TRUE)[,1]
  frac_invalid <- ifelse(count_forest > 0, count_true / count_forest, NA_real_)

  stats_df <- data.frame(
    year        = yrs_lab,
    mean_sum    = mean_sum,
    frac_invalid_gt_tol = frac_invalid
  )

  # Write results
  write.csv(stats_df,
            file.path(out_dir, "Summary_Stats_forestOnly.csv"),
            row.names = FALSE)

  # Quick console report
  message(glue("Years checked: {min(yrs_lab)}–{max(yrs_lab)} ({length(yrs_lab)} years)"))
  message(glue("Mean(mean_sum) inside forest: {round(mean(stats_df$mean_sum, na.rm=TRUE), 6)}"))
  message(glue("Max |mean_sum-1| inside forest: {round(max(abs(stats_df$mean_sum-1), na.rm=TRUE), 6)}"))
  message(glue("Median frac invalid (> {tolerance}) inside forest: ",
               "{scales::percent(median(stats_df$frac_invalid_gt_tol, na.rm=TRUE))}"))

  # Optional diagnostic plot saved to disk
  p1 <- ggplot(stats_df, aes(year, mean_sum)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(title = glue("Forest-only Mean of (aff+mgd) — {scenario}"),
         y = "Mean (aff+mgd) inside forest") +
    theme_bw()

  p2 <- ggplot(stats_df, aes(year, frac_invalid_gt_tol)) +
    geom_line() +
    labs(title = glue("Share of forest cells with |aff+mgd-1| > {tolerance}"),
         y = "Fraction invalid") +
    theme_bw()

  ggsave(file.path(out_dir, glue("Sanity_{scenario}_meanSum_forestOnly.png")), p1,
         width = 10, height = 4, dpi = 220)
  ggsave(file.path(out_dir, glue("Sanity_{scenario}_invalidFrac_forestOnly.png")), p2,
         width = 10, height = 4, dpi = 220)

  message("Sanity check outputs written to: ", out_dir)
}

message("\nAll done.")
