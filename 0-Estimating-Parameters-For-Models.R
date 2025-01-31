# Estimating parameters

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, sn, moments)

# kgc loads depedencies that conflict with tidyverse call kgc:: instead

# Import data ------------------------------------------------------------------
data <- read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS_20240807.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)


gauge_information <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS_20240807.csv",
  show_col_types = FALSE
) |> 
  select(gauge, state, lat, lon)



# Assign climate type to each gauge using `kgc` --------------------------------
## `LookupCZ` function provides the climate zone based on lon and lat 
## Relies on the climatezone dataframe
climatezones <- kgc::climatezones

## To use LookupCZ the data must be in |site_ID|lon|lat| format ================
formatted_gauge_information <- gauge_information |> 
  select(!state) |> 
  relocate(
    lon,
    .after = 1
  ) |> 
  mutate(
    rndCoord.lon = kgc::RoundCoordinates(lon),
    rndCoord.lat = kgc::RoundCoordinates(lat)
  ) 


## Gauge information with climate type =========================================
climate_type_gauge_info <- cbind(
  formatted_gauge_information, 
  "climate_type" = kgc::LookupCZ(data = formatted_gauge_information)
) |> 
  as_tibble() |> 
  mutate(
    major_climate_type = str_sub(climate_type, start = 1L, end = 1L)
  ) |> 
  select(gauge, major_climate_type) |> 
  # Nice names for plotting
  mutate(
    major_climate_type = case_when(
      major_climate_type == "A" ~ "Tropical (A)",
      major_climate_type == "B" ~ "Dry (B)",
      major_climate_type == "C" ~ "Temperate (C)",
      .default = major_climate_type
    )
  ) |>
  # Make climate types in set order
  mutate(
    major_climate_type = factor(
      major_climate_type, 
      levels = c("Overall", "Tropical (A)", "Dry (B)", "Temperate (C)")
      )
  )




# Rainfall mean, sd, auto and skewness -----------------------------------------
get_lag_1_autocorrelation <- function(timeseries) {
  acf(timeseries, na.action = na.pass, plot = FALSE)$acf[2]
}

rainfall_stats_per_gauge <- data |>
  summarise(
    p_mean = mean(p_mm),
    p_sd = sd(p_mm),
    p_auto = get_lag_1_autocorrelation(p_mm),
    p_skew = skewness(p_mm),
    .by = gauge
  ) |> 
  left_join(
    climate_type_gauge_info,
    by = join_by(gauge)
  )


## Overall rainfall statistics =================================================
summary_rainfall_stat <- rainfall_stats_per_gauge |>
  pivot_longer(
    cols = starts_with("p"),
    names_to = "metric",
    values_to = "values"
  ) 

overall_summary_rainfall_stat <- summary_rainfall_stat |>
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  )


## By climate type rainfall statistics =========================================
climate_type_summary_rainfall_stat <- summary_rainfall_stat |> 
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = c(metric, major_climate_type)
  )


## Boxplot of rainfall statistics ==============================================
### To add a overall to boxplot copy the entire tibble replace major climate 
### type with overall and rbind() (not very elegant)
plot_summary_rainfall_stat <- summary_rainfall_stat |> 
  mutate(
    major_climate_type = as.factor("Overall") # keep major_climate_type as factor
  ) |> 
  rbind(
    summary_rainfall_stat
  ) |> 
  # give metric a nice name
  mutate(
    metric = case_when(
      metric == "p_auto" ~ "Autocorrelation",
      metric == "p_mean" ~ "Mean",
      metric == "p_sd" ~ "Standard Deviation",
      metric == "p_skew" ~ "Skewness",
      .default = metric
    )
  ) |> 
  # order metric
  mutate(
    metric = factor(
      metric, 
      levels = c("Mean", "Standard Deviation", "Autocorrelation", "Skewness")
      )
  )



### Code to produce n = ...
###  Not great it needs the have mulitiple y_positions for each facet
### Another alternative is to have it only for the first facet?
### FIX LATER
count_summary_rainfall_data <- plot_summary_rainfall_stat |> 
  summarise(
    y_position = 0,
    n = n(),
    .by = major_climate_type
  ) |> 
  mutate(
    label_n = paste0("n = ", n)
  )


### Plot #######################################################################
plot_summary_rainfall_stat |> 
  ggplot(aes(x = major_climate_type, y = values, fill = major_climate_type)) +
  geom_boxplot(show.legend = FALSE) +
  geom_text(
    aes(x = major_climate_type, y = y_position, label = label_n),
    data = count_summary_rainfall_data,
    inherit.aes = FALSE
    ) +
  labs(
    x = "Major Climate Type",
    y = "Value"
  ) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme()


# Calculate average fitted slope and intercept ---------------------------------
get_fitted_intercept <- function(rainfall, streamflow) {
  coef(lm(streamflow ~ rainfall))[1]
}

get_fitted_slope <- function(rainfall, streamflow) {
  coef(lm(streamflow ~ rainfall))[2]
}

line_of_best_fit_per_catchment <- data |>
  summarise(
    fitted_intercept = get_fitted_intercept(p_mm, bc_q),
    fitted_slope = get_fitted_slope(p_mm, bc_q),
    .by = gauge
  ) 


## Overall intercept and slope statistics ======================================
summary_intercept_and_slope <- line_of_best_fit_per_catchment |>
  pivot_longer(
    cols = starts_with("fitted"),
    names_to = "metric",
    values_to = "values"
  ) 

overall_summary_intercept_and_slope <- summary_intercept_and_slope |>
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  )


## By climate type intercept and slope statistics ==============================
### ADD HERE

# Calculate average autocorrelation --------------------------------------------
summary_autocorrelation <- data |> 
  summarise(
    auto = get_lag_1_autocorrelation(bc_q),
    .by = gauge
    ) |> 
  # Make it the same format as other code
  pivot_longer(
    cols = auto,
    names_to = "metric",
    values_to = "values"
  )


## Overall autocorrelation statistics ==========================================
overall_summary_autocorrelation <- summary_autocorrelation |> 
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  ) 


## By climate type autocorrelation statistics ==================================


# Calculate average sd and skewness of residuals around line of best fit -------
get_sd_around_line_of_best_fit <- function(rainfall, streamflow) {
  sd(lm(streamflow ~ rainfall)$residuals)
}

get_skew_around_line_of_best_fit <- function(rainfall, streamflow) {
  skewness(lm(streamflow ~ rainfall)$residuals)
}

spread_and_shape_around_line_of_best_fit <- data |>
  summarise(
    sd_residuals = get_sd_around_line_of_best_fit(p_mm, bc_q),
    skew_residuals = get_skew_around_line_of_best_fit(p_mm, bc_q),
    .by = gauge
  )


## Overall spread and shape statistics =========================================

summary_spread_and_shape <- spread_and_shape_around_line_of_best_fit |>
  pivot_longer(
    cols = ends_with("residuals"),
    names_to = "metric",
    values_to = "values"
  ) 

overall_summary_spread_and_shape <- summary_spread_and_shape |> 
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  )


## By climate type spread and shape statistics =================================

# Combine streamflow model parameters into a single table ----------------------
## Overall =====================================================================
summary_streamflow_model <- rbind(
  summary_intercept_and_slope,
  summary_autocorrelation,
  summary_spread_and_shape
)


## By climate type =============================================================
### combine everything for plotting similar table to the rainfall one
summary_partitioning_stat <- rbind(
  summary_intercept_and_slope,
  summary_autocorrelation,
  summary_spread_and_shape
) |> 
  left_join(
    climate_type_gauge_info,
    by = join_by(gauge)
  )

### Again, this is to add Overall to facets
plot_summary_partitioning_stat <- summary_partitioning_stat |> 
  mutate(
    major_climate_type = as.factor("Overall")
  ) |> 
  rbind(
    summary_partitioning_stat
  ) |> 
  # Rename facets
  mutate(
    metric = case_when(
      metric == "auto" ~ "Autocorrelation",
      metric == "fitted_intercept" ~ "Intercept",
      metric == "fitted_slope" ~ "Slope",
      metric == "sd_residuals" ~ "Standard Deviation",
      metric == "skew_residuals" ~ "Skewness",
      .default = metric
    )
  ) |> 
  # Order facets
  mutate(
    metric = factor(
      metric, 
      levels = c("Intercept", "Slope", "Autocorrelation", "Standard Deviation", "Skewness")
      )
  )


### Plot #######################################################################
plot_summary_partitioning_stat |> 
  ggplot(aes(x = major_climate_type, y = values, fill = major_climate_type)) +
  geom_boxplot(show.legend = FALSE) +
  labs(
    x = "Major Climate Type",
    y = "Value"
  ) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme()

