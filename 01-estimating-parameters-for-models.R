# Estimating parameters

# TODO:
# 3. Fix up n = count for the boxplots
# 4. Try scatter plot method




# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")


# Import libraries -------------------------------------------------------------
library(tidyverse)
library(sn)
library(moments)


# 1. Import data ------------------------------------------------------------------
data <- read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)


gauge_info <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
) |>
  select(!c(station_name, bc_lambda)) |> 
  # Make climate types in set order
  mutate(
    major_climate_type = factor(
      major_climate_type, 
      levels = c("Overall", "Tropical (A)", "Dry (B)", "Temperate (C)")
    )
  )


# 2. Examine variety of climate types ------------------------------------------
number_climate_types <- gauge_info |> 
  summarise(
    n = n(),
    .by = major_climate_type
  )

# Temperate = 476, Dry = 48 and Tropical = 37
# Focus on temperate catchments


# 3. Rainfall statistics -------------------------------------------------------
## 3.1 mean, sd, auto and skewness calculation =================================
get_lag_1_autocorrelation <- function(timeseries) {
  acf(timeseries, na.action = na.pass, plot = FALSE)$acf[2]
}

rainfall_stats_per_gauge <- data |>
  summarise(
    "Mean" = mean(p_mm),
    "Standard Deviation" = sd(p_mm),
    "Autocorrelation" = get_lag_1_autocorrelation(p_mm),
    "Skewness" = skewness(p_mm),
    .by = gauge
  ) |> 
  left_join(
    gauge_info,
    by = join_by(gauge)
  )


summary_rainfall_stat <- rainfall_stats_per_gauge |>
  pivot_longer(
    cols = Mean:Skewness,
    names_to = "metric",
    values_to = "values"
  ) 


## 3.2 Determine rainfall control and multiplier statistics =====================

rainfall_control_and_multipliers <- summary_rainfall_stat |> 
  summarise(
    median = median(values),
    q1 = quantile(values, 0.01, names = FALSE),
    q25 = quantile(values, 0.25, names = FALSE),
    q75 = quantile(values, 0.75, names = FALSE),
    q99 = quantile(values, 0.99, names = FALSE),
    .by = c(metric, major_climate_type)
  ) |> 
  mutate(
    q1_multi = q1 / median,
    q25_multi = q25 / median,
    q75_multi = q75 / median,
    q99_multi = q99 / median
  )


## Direction of changes for multiplier selection (from literature):
### - Mean = decrease
### - Standard deviation = increase
### - Autocorrelation = increase
### - Skewness = increase

rainfall_temperate_parameters <- rainfall_control_and_multipliers |> 
  filter(major_climate_type == "Temperate (C)") |> 
  filter(metric != "Mean") |> 
  select(metric, major_climate_type, median, q75, q99, q75_multi, q99_multi) |> 
  rename(
    control_parameter = median,
    small_change_parameter = q75,
    large_change_parameter = q99,
    small_multiplier = q75_multi,
    large_multiplier = q99_multi
  )

### Mannually fix the mean rainfall so its a decrease
just_mean_rainfall <- rainfall_control_and_multipliers |> 
  filter(major_climate_type == "Temperate (C)") |> 
  filter(metric == "Mean") |> 
  select(!c(q75, q99, q75_multi, q99_multi)) |> 
  rename(
    control_parameter = median,
    small_change_parameter = q25,
    large_change_parameter = q1,
    small_multiplier = q25_multi,
    large_multiplier = q1_multi
  )

rainfall_temperate_parameters <- rbind(just_mean_rainfall, rainfall_temperate_parameters) |> 
  # order metric
  mutate(
    metric = factor(
      metric, 
      levels = c("Mean", "Standard Deviation", "Autocorrelation", "Skewness")
    )
  )


write_csv(
  rainfall_temperate_parameters,
  "./Results/rainfall_temperate_parameters.csv"
  )

## 3.3 Tidy summary stats ready for plotting ===================================
### To add a overall to boxplot copy the entire tibble replace major climate 
### type with overall and rbind() (not very elegant)
plot_summary_rainfall_stat <- summary_rainfall_stat |> 
  mutate(
    major_climate_type = as.factor("Overall") # keep major_climate_type as factor
  ) |> 
  rbind(
    summary_rainfall_stat
  ) |> 
  # order metric
  mutate(
    metric = factor(
      metric, 
      levels = c("Mean", "Standard Deviation", "Autocorrelation", "Skewness")
      )
  )




## 3.4 Boxplot of rainfall statistics ==========================================
rainfall_boxplot <- plot_summary_rainfall_stat |> 
  ggplot(aes(x = major_climate_type, y = values, fill = major_climate_type)) +
  geom_boxplot(
    show.legend = FALSE,
    outliers = FALSE,
    staplewidth = 0.25
    ) +
  geom_hline(
    aes(yintercept = small_change_parameter), 
    data = rainfall_temperate_parameters |> filter(major_climate_type == "Temperate (C)"),
    linetype = "dashed",
    colour = "#ff7f00",
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = large_change_parameter), 
    data = rainfall_temperate_parameters |> filter(major_climate_type == "Temperate (C)"),
    linetype = "dotdash",
    colour = "#a65628",
    linewidth = 1
  ) +
  labs(
    x = "Major Climate Type",
    y = "Value"
  ) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme()


### Save graph #################################################################
ggsave(
  filename = "rainfall_boxplot.pdf",
  plot = rainfall_boxplot,
  device = "pdf",
  path = "./Graphs",
  width = 297,
  height = 210,
  units = "mm"
)



# 4. Rainfall-partitioning statistics ------------------------------------------

## 4.1 Calculate fitted slope and intercept per gauge ==========================
get_fitted_intercept <- function(rainfall, streamflow) {
  coef(lm(streamflow ~ rainfall))[1]
}

get_fitted_slope <- function(rainfall, streamflow) {
  coef(lm(streamflow ~ rainfall))[2]
}

line_of_best_fit_per_catchment <- data |>
  summarise(
    Intercept = get_fitted_intercept(p_mm, bc_q),
    Slope = get_fitted_slope(p_mm, bc_q),
    .by = gauge
  ) 


summary_intercept_and_slope <- line_of_best_fit_per_catchment |>
  pivot_longer(
    cols = !gauge,
    names_to = "metric",
    values_to = "values"
  ) 



## 4.2 Calculate boxcox streamflow autocorrelation per gauge ===================
summary_autocorrelation <- data |> 
  summarise(
    Autocorrelation = get_lag_1_autocorrelation(bc_q),
    .by = gauge
    ) |> 
  # Make it the same format as other code
  pivot_longer(
    cols = Autocorrelation,
    names_to = "metric",
    values_to = "values"
  )



## 4.3 Calculate sd and skewness of residuals around line of best fit per gauge ====
get_sd_around_line_of_best_fit <- function(rainfall, streamflow) {
  sd(lm(streamflow ~ rainfall)$residuals)
}

get_skew_around_line_of_best_fit <- function(rainfall, streamflow) {
  skewness(lm(streamflow ~ rainfall)$residuals)
}

spread_and_shape_around_line_of_best_fit <- data |>
  summarise(
    "Standard Deviation" = get_sd_around_line_of_best_fit(p_mm, bc_q),
    "Skewness" = get_skew_around_line_of_best_fit(p_mm, bc_q),
    .by = gauge
  )


summary_spread_and_shape <- spread_and_shape_around_line_of_best_fit |>
  pivot_longer(
    cols = !gauge,
    names_to = "metric",
    values_to = "values"
  ) 


## 4.4 Bring intercept, slope, auto, sd and skew results together ==============

summary_partitioning <- rbind(
  summary_intercept_and_slope,
  summary_autocorrelation,
  summary_spread_and_shape
  ) |>
left_join(
  gauge_info,
  by = join_by(gauge)
  )



## 4.5 Determine partitioning control and multiplier statistics ================

partioning_control_and_multipliers <- summary_partitioning |> 
  summarise(
    median = median(values),
    q1 = quantile(values, 0.01, names = FALSE),
    q25 = quantile(values, 0.25, names = FALSE),
    q75 = quantile(values, 0.75, names = FALSE),
    q99 = quantile(values, 0.99, names = FALSE),
    .by = c(metric, major_climate_type)
  ) |> 
  mutate(
    q1_multi = q1 / median,
    q25_multi = q25 / median,
    q75_multi = q75 / median,
    q99_multi = q99 / median
  )


## Direction of changes for multiplier selection (from literature):
### - Intercept = decrease
### - Slope = increase 
### - Autocorrelation = increase
### - Standard deviation = increase
### - Skewness = increase

partitioning_temperate_parameters <- partioning_control_and_multipliers |> 
  filter(major_climate_type == "Temperate (C)") |> 
  filter(metric != "fitted_intercept") |> 
  select(metric, major_climate_type, median, q75, q99, q75_multi, q99_multi) |> 
  rename(
    control_parameter = median,
    small_change_parameter = q75,
    large_change_parameter = q99,
    small_multiplier = q75_multi,
    large_multiplier = q99_multi
  )

### Mannually fix the mean rainfall so its a decrease
just_intercept <- partioning_control_and_multipliers |> 
  filter(major_climate_type == "Temperate (C)") |> 
  filter(metric == "fitted_intercept") |> 
  select(!c(q75, q99, q75_multi, q99_multi)) |> 
  rename(
    control_parameter = median,
    small_change_parameter = q25,
    large_change_parameter = q1,
    small_multiplier = q25_multi,
    large_multiplier = q1_multi
  )

partitioning_temperate_parameters <- rbind(just_intercept, partitioning_temperate_parameters) |> 
  # Order facets
  mutate(
    metric = factor(
      metric, 
      levels = c("Intercept", "Slope", "Autocorrelation", "Standard Deviation", "Skewness")
    )
  )



write_csv(
  partitioning_temperate_parameters,
  "./Results/partitioning_temperate_parameters.csv"
)


### Again, this is to add Overall to facets
plot_summary_partitioning_stat <- summary_partitioning |> 
  mutate(
    major_climate_type = as.factor("Overall")
  ) |> 
  rbind(
    summary_partitioning
  ) |> 
  # Order facets
  mutate(
    metric = factor(
      metric, 
      levels = c("Intercept", "Slope", "Autocorrelation", "Standard Deviation", "Skewness")
      )
  )



## 4.6 Boxplot of rainfall-partitioning statistics =============================
partitioning_boxplot <- plot_summary_partitioning_stat |> 
  ggplot(aes(x = major_climate_type, y = values, fill = major_climate_type)) +
  geom_boxplot(
    show.legend = FALSE,
    outliers = FALSE,
    staplewidth = 0.25
    ) +
  geom_hline(
    aes(yintercept = small_change_parameter), 
    data = partitioning_temperate_parameters |> filter(major_climate_type == "Temperate (C)"),
    linetype = "dashed",
    colour = "#ff7f00",
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = large_change_parameter), 
    data = partitioning_temperate_parameters |> filter(major_climate_type == "Temperate (C)"),
    linetype = "dotdash",
    colour = "#a65628",
    linewidth = 1
  ) +
  labs(
    x = "Major Climate Type",
    y = "Value"
  ) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme()


### Save plot ##################################################################
ggsave(
  filename = "partitioning_boxplot.pdf",
  plot = partitioning_boxplot,
  device = "pdf",
  path = "./Graphs",
  width = 297,
  height = 210,
  units = "mm"
)
