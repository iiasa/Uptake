library(tidyverse)
library(lubridate)
library(ggthemes)
library(scales)

# Define file paths
base_path <- "H:/UPTAKE/burned_area_csv/"
scenarios <- c("historical", "ssp126", "ssp245", "ssp370")

# File lists
gfdl_files <- paste0(base_path, c(
  "gfdl_historical.csv",
  "gfdl_ssp126.csv",
  "gfdl_ssp245.csv",
  "gfdl_ssp370.csv"
))

ukesm_files <- paste0(base_path, c(
  "ukesm_Historical.csv",
  "ukesm_SSP126.csv",
  "ukesm_SSP245.csv",
  "ukesm_SSP370.csv"
))

# Robust data processing function
process_burned_data <- function(file, scenario, model) {
  # Read the entire file as text
  raw_text <- readLines(file, warn = FALSE)
  
  # Extract the three main lines
  time_line <- raw_text[1]
  obs_line <- raw_text[2]
  pred_line <- raw_text[3]
  
  # Split lines into values
  time_vals <- str_split(time_line, ",")[[1]] %>% 
    str_trim() %>% 
    .[-1]  # Remove the first element which is likely a header
  
  obs_vals <- str_split(obs_line, ",")[[1]] %>% 
    str_trim() %>% 
    .[-1]  # Remove the first element which is likely a header
  
  pred_vals <- str_split(pred_line, ",")[[1]] %>% 
    str_trim() %>% 
    .[-1]  # Remove the first element which is likely a header
  
  # Create data frame
  monthly_data <- tibble(
    time_str = time_vals,
    Observation = as.numeric(obs_vals),
    Prediction = as.numeric(pred_vals),
    scenario = scenario,
    model = model
  ) %>%
    filter(nchar(time_str) == 6) %>%  # Keep only monthly data (YYYYMM format)
    mutate(
      date = ymd(paste0(time_str, "01")),
      year = year(date),
      month = month(date)
    )
  
  # Create yearly data
  yearly_data <- monthly_data %>%
    group_by(year, scenario, model) %>%
    summarize(
      Observation_sum = sum(Observation, na.rm = TRUE),
      Prediction_sum = sum(Prediction, na.rm = TRUE),
      Prediction_sd = sd(Prediction, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(yearly_data)
}

# Safe processing function
process_files_safely <- function(files, scenarios, model_name) {
  map2_dfr(files, scenarios, ~{
    tryCatch({
      process_burned_data(.x, .y, model_name)
    }, error = function(e) {
      message(paste("Error processing", .x, ":", e$message))
      return(NULL)
    })
  })
}

# Process GFDL files
gfdl_data <- process_files_safely(gfdl_files, scenarios, "GFDL")

# Process UKESM files
ukesm_data <- process_files_safely(ukesm_files, scenarios, "UKESM")

# Combine data
all_data <- bind_rows(gfdl_data, ukesm_data)

# Prepare data for plotting - connect historical to ssp126
plot_data <- all_data %>%
  # Create continuous series for each scenario
  mutate(
    continuous_scenario = case_when(
      scenario == "historical" ~ "Historical",
      scenario == "ssp126" ~ "Mit2p6",
      scenario == "ssp245" ~ "Mit4p5",
      scenario == "ssp370" ~ "Mi7p0"
    )
  ) %>%
  # Create a special connection between historical and ssp126
  mutate(
    connected_scenario = case_when(
      scenario == "historical" ~ "Historical-Mit2p6",
      scenario == "ssp126" ~ "Historical-Mit2p6",
      TRUE ~ continuous_scenario
    )
  ) %>%
  filter(model == "UKESM") %>%
  select(year, scenario, continuous_scenario, connected_scenario, Prediction_sum) %>%
  left_join(
    all_data %>%
      filter(model == "GFDL") %>%
      mutate(
        continuous_scenario = case_when(
          scenario == "historical" ~ "Historical",
          scenario == "ssp126" ~ "Mit2p6",
          scenario == "ssp245" ~ "Mit4p5",
          scenario == "ssp370" ~ "Mi7p0"
        ),
        connected_scenario = case_when(
          scenario == "historical" ~ "Historical-Mit2p6",
          scenario == "ssp126" ~ "Historical-Mit2p6",
          TRUE ~ continuous_scenario
        )
      ) %>%
      select(year, scenario, sd = Prediction_sd),
    by = c("year", "scenario")
  ) %>%
  mutate(
    ymin = pmax(Prediction_sum - sd, 0.1),
    ymax = Prediction_sum + sd
  )

# Get observation data from ALL scenarios
obs_data <- all_data %>%
  select(year, scenario, Observation_sum) %>%
  distinct() %>%
  filter(!is.na(Observation_sum), Observation_sum > 0) %>%
  mutate(scenario = case_when(
    scenario == "historical" ~ "Historical",
    scenario == "ssp126" ~ "Mit2p6",
    scenario == "ssp245" ~ "Mit4p5",
    scenario == "ssp370" ~ "Mi7p0"
  ))

# Prepare line data with separated segments
line_data <- bind_rows(
  # Historical segment (blue)
  plot_data %>% 
    filter(continuous_scenario == "Historical", year <= 2014) %>%
    mutate(line_color = "Historical"),
  # Mit2p6 segment (green)
  plot_data %>% 
    filter(continuous_scenario == "Mit2p6", year >= 2015) %>%
    mutate(line_color = "Mit2p6"),
  # Other scenarios
  plot_data %>% 
    filter(!continuous_scenario %in% c("Historical", "Mit2p6")) %>%
    mutate(line_color = continuous_scenario)
)

# Create color palette
scenario_colors_line <- c(
  "Historical" = "#1f77b4",       # Blue for historical
  "Mit2p6" = "#2ca02c",           # Green for Mit2p6
  "Mit4p5" = "#ff7f0e",           # Orange
  "Mi7p0" = "#d62728",            # Red
  "Observations" = "black"        # Black for observations
)

scenario_colors_fill <- scenario_colors_line[1:4]

# Create plot with integrated legend
ggplot() +
  # SD Ribbon
  geom_ribbon(
    data = plot_data,
    aes(x = year, ymin = ymin, ymax = ymax, fill = continuous_scenario),
    alpha = 0.3
  ) +
  # Main UKESM lines - with separated segments
  geom_line(
    data = line_data,
    aes(x = year, y = Prediction_sum, color = line_color, group = connected_scenario),
    size = 1.2
  ) +
  # Observation data from ALL scenarios
  geom_line(
    data = obs_data,
    aes(x = year, y = Observation_sum, color = "Observations"),
    size = 1.2
  ) +
  geom_point(
    data = obs_data,
    aes(x = year, y = Observation_sum, color = "Observations"),
    size = 2.5,
    shape = 21,
    fill = "white"
  ) +
  # Vertical line at 2015
  geom_vline(xintercept = 2015, linetype = "dashed", color = "gray30", size = 0.8) +
  # Unified color and fill scales
  scale_color_manual(
    name = "Scenario",
    values = scenario_colors_line,
    breaks = c("Historical", "Mit2p6", "Mit4p5", "Mi7p0", "Observations"),
    labels = c("Historical", "Mit2p6", "Mit4p5", "Mi7p0", "Observations")
  ) +
  scale_fill_manual(
    name = "Scenario",
    values = scenario_colors_fill,
    breaks = c("Historical", "Mit2p6", "Mit4p5", "Mi7p0"),
    labels = c("Historical", "Mit2p6", "Mit4p5", "Mi7p0")
  ) +
  labs(
    title = "Burned Area Projections",
    x = "Year",
    y = "Burned Area (Annual Sum)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c(rep(1, 4), 1),
        shape = c(rep(NA, 4), 21),
        fill = c(rep(NA, 4), "white")
      )
    ),
    fill = guide_legend(override.aes = list(alpha = 0.3))
  ) +
  # Historical period annotation
  annotate("text", x = 2000, y = max(plot_data$ymax) * 0.95, 
           label = "Historical Period", color = "gray30", size = 4) +
  # Future period annotation
  annotate("text", x = 2050, y = max(plot_data$ymax) * 0.95, 
           label = "Future Projections", color = "gray30", size = 4)

# Save plot
ggsave("burned_area_final.png", width = 10, height = 7, dpi = 300)