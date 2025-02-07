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
    # this probably should be map then expand the list 
    # rather than copy and pasting the same function many times
    median = median(values),
    q5 = quantile(values, 0.05, names = FALSE), 
    q25 = quantile(values, 0.25, names = FALSE),
    q75 = quantile(values, 0.75, names = FALSE),
    q95 = quantile(values, 0.95, names = FALSE),
    q99 = quantile(values, 0.99, names = FALSE),
    .by = c(metric, major_climate_type)
  ) 


### make_control_and_multiplier_parameters function ############################
# example inputs:
# metric = "Intercept"
# small_change = "q5"
# large_change = "q1"

make_control_and_multiplier_parameters <- function(metric, small_change, large_change, data) {
  data |> 
    filter(metric == {{ metric }}) |> 
    select(median, {{ small_change }}, {{ large_change }}) |>
    rename(
      small_change = {{ small_change }},
      large_change = {{ large_change }}
    ) |> 
    mutate(
      small_multiplier = small_change / median,
      large_multiplier = large_change / median
    ) |> 
    add_column(
      metric = {{ metric }},
      .before = 1
    )
}


## Direction of changes for multiplier selection (from literature):
### - Mean = decrease
### - Standard deviation = increase
### - Autocorrelation = increase
### - Skewness = increase


temperate_rainfall_control_and_multipliers <- rainfall_control_and_multipliers |> 
  filter(major_climate_type == "Temperate (C)")

rainfall_temperate_parameters <- pmap(
  .l = list(
    temperate_rainfall_control_and_multipliers |> pull(metric),
    # Mean small q25, large q5
    # SD small q75, large q95
    # Autocorrelation small q75, large q95
    # Skewness small q75, large q95
    c("q25", "q75", "q75", "q75"), # small change
    c("q5", "q95", "q99", "q99") # large change
  ),
  .f = make_control_and_multiplier_parameters,
  data = temperate_rainfall_control_and_multipliers
) |> 
  list_rbind() |> 
  rename(control = median) |> 
  # Order facets
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
    aes(yintercept = small_change), 
    data = rainfall_temperate_parameters,
    linetype = "dashed",
    colour = "#ff7f00",
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = large_change), 
    data = rainfall_temperate_parameters,
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
### Progressively step up percentiles until we see a change:
### - Start at 25th/75th (depending or direction we want),
### - Step up to 80th (sd big change)
### - Step up to 90th
### - Step up to 95th (skew big change)
### - Step up to 99th (auto-big change can't be cause skewness exceeds limit)



partioning_control_and_multipliers <- summary_partitioning |> 
  summarise(
    # this probably should be map then expand the list 
    # rather than copy and pasting the same function many times
    median = median(values),
    q5 = quantile(values, 0.05, names = FALSE), 
    q25 = quantile(values, 0.25, names = FALSE),
    q65 = quantile(values, 0.65, names = FALSE),
    q75 = quantile(values, 0.75, names = FALSE),
    q80 = quantile(values, 0.8, names = FALSE),
    q95 = quantile(values, 0.95, names = FALSE),
    q99 = quantile(values, 0.99, names = FALSE),
    .by = c(metric, major_climate_type)
  ) 


## Direction of changes for multiplier selection (from literature):
### - Intercept = decrease
### - Slope = increase 
### - Autocorrelation = increase
### - Standard deviation = increase
### - Skewness = increase


temperate_partioning_control_and_multipliers <- partioning_control_and_multipliers |> 
  filter(major_climate_type == "Temperate (C)")

partitioning_temperate_parameters <- pmap(
  .l = list(
    temperate_partioning_control_and_multipliers |> pull(metric),
    # Interecpt small q25
    # Slope small q75
    # Auto large q99
    # Sd large q80
    # Skew large q95 (can't do q99 due to skewness limits)
    c("q25", "q65", "q75", "q75", "q75"), # small change
    c("q5", "q95", "q99", "q80", "q95") # large change
  ),
  .f = make_control_and_multiplier_parameters,
  data = temperate_partioning_control_and_multipliers
) |> 
  list_rbind() |> 
  rename(control = median) |> 
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

### Is there a catchment with values similar to medians? #######################
median_partitioning_parameters <- partitioning_temperate_parameters |> 
  select(metric, control) |> 
  mutate(
    metric = as.character(metric)
  )


find_catchments_near_medians <- function(metric, metric_median, tol, data) {
  data |> 
    filter(metric == {{ metric }}) |> 
    filter(near(values, as.numeric(metric_median), tol)) |> 
    pull(gauge)
}


gauges_near_medians <- pmap(
  .l = list(
    median_partitioning_parameters |> pull(metric),
    median_partitioning_parameters |> pull(control),
    c(0.25, 0.01, 0.1, 0.5, 0.1) # play with tols to get smallest number catchments
  ),
  .f = find_catchments_near_medians,
  data = summary_partitioning
) |> 
  reduce(
    .f = base::intersect # find common elements in vectors
  )

representative_catchments <- summary_partitioning |> 
  filter(gauge %in% gauges_near_medians) |> 
  arrange(gauge)






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
    aes(yintercept = small_change), 
    data = partitioning_temperate_parameters,
    linetype = "dashed",
    colour = "#ff7f00",
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = large_change), 
    data = partitioning_temperate_parameters,
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

