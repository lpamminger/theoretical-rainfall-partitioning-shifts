# Estimating parameters

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")




## PUT CLIMATE ANALYSIS HERE ##
# method:
# 1. load kgc package
# 2. Get the gauge_information tibble into |siteid|lon|lat
# 3. Round lat and lon using RoundCoordinates (use specified names)
# 4. Use LookupCZ to turn lat+lon into climate type
# 5. Save result
# 6. Repeat analysis below by climate type


## TEMP CODE transferred from testing ## - REMOVE

library(kgc)
library(tidyverse)

round_any = function(x, accuracy, f = round){
  f(x / accuracy) * accuracy
}

gauge_information <- read_csv(
  "PHD/Papers/RQ1/theoretical-rainfall-partitioning-shifts/Data/Tidy/gauge_information_CAMELS_20240807.csv",
  show_col_types = FALSE
) |> 
  select(gauge, state, lat, lon)

climate_zones <- climatezones |> 
  as_tibble() |> 
  rename(
    lat = Lat,
    lon = Lon
  )



data <- data.frame(Site = c("GC","UFS","NEG"),
                   Longitude = c(-15.42,10.98,34.78),
                   Latitude = c(27.82,47.42,30.86))
data <- data.frame(data,
                   rndCoord.lon = RoundCoordinates(data$Longitude),
                   rndCoord.lat = RoundCoordinates(data$Latitude))
data <- data.frame(data,ClimateZ=LookupCZ(data))





# use the function?
# data must be |site_ID|lon|lat|
adjusted_gauge_information <- gauge_information |> 
  select(!state) |> 
  relocate(
    lon,
    .after = 1
  ) |> 
  mutate(
    rndCoord.lon = RoundCoordinates(lon),
    rndCoord.lat = RoundCoordinates(lat)
  ) #|> 
#as.data.frame()




x <- cbind(
  adjusted_gauge_information, 
  LookupCZ(data = adjusted_gauge_information)
) |> 
  as_tibble() |> 
  rename(
    climate_type = `LookupCZ(data = adjusted_gauge_information)`
  ) |> 
  mutate(
    major_climate_type = str_sub(climate_type, start = 1L, end = 1L)
  )







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
