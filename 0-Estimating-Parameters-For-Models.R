# Estimating parameters

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, sn, moments)

# Import data ------------------------------------------------------------------
data <- read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS_20240807.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
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
  ) 

summary_rainfall_stat <- rainfall_stats_per_gauge |>
  pivot_longer(
    cols = starts_with("p"),
    names_to = "metric",
    values_to = "values"
  ) |>
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  )


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

summary_intercept_and_slope <- line_of_best_fit_per_catchment |>
  pivot_longer(
    cols = starts_with("fitted"),
    names_to = "metric",
    values_to = "values"
  ) |>
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  )


# Calculate average autocorrelation --------------------------------------------
autocorrelation_per_catchment <- data |> 
                                   summarise(
                                     auto = get_lag_1_autocorrelation(bc_q),
                                     .by = gauge
                                   ) 

summary_autocorrelation <- autocorrelation_per_catchment |> 
  summarise(
    q50 = median(auto),
    q5 = quantile(auto, 0.05),
    q95 = quantile(auto, 0.95),
  ) |> 
  add_column(
    "metric" = "auto",
    .before = 1
  )



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

summary_spread_and_shape <- spread_and_shape_around_line_of_best_fit |>
  pivot_longer(
    cols = ends_with("residuals"),
    names_to = "metric",
    values_to = "values"
  ) |> 
  summarise(
    q50 = median(values),
    q5 = quantile(values, 0.05),
    q95 = quantile(values, 0.95),
    .by = metric
  )

# Combine streamflow model parameters into a single table ----------------------
summary_streamflow_model <- rbind(
  summary_intercept_and_slope,
  summary_autocorrelation,
  summary_spread_and_shape
)
