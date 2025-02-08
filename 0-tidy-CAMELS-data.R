# Preparing CAMEL-AUS v2 data

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")


# Import libraries -------------------------------------------------------------
library(MASS) # boxcox conversion
library(tidyverse)
library(checkmate) # start stop indexes
# kgc used for finding climate zones 
# kgc loads dependencies that conflict with tidyverse call kgc:: instead

# Functions for tidying --------------------------------------------------------
source("./Functions/boxcox_transforms.R")
source("./Functions/utility.R")

# TODO: # remove gauges with less least 10 years of continuous data

# 1. Import CAMELS v2 data -----------------------------------------------------

daily_streamflow_raw <- readr::read_csv(
  "./Data/Raw/streamflow_mmd.csv",
  na = c("-99.99"), # NA is represented using -99.99 in CAMELS dataset. Convert to NA
  show_col_types = FALSE # quiets column specification
)

daily_precip_raw <- readr::read_csv(
  "./Data/Raw/precipitation_AGCD.csv",
  na = c("-99.99"),
  show_col_types = FALSE
)

catchment_information <- readr::read_csv(
  "./Data/Raw/CAMELS_AUS_Attributes&Indices_MasterTable.csv",
  col_select = c(
    "station_id",
    "station_name",
    "state_outlet",
    "lat_outlet",
    "long_outlet"
  ),
  show_col_types = FALSE
) |> 
  rename(
    gauge = station_id,
    lat = lat_outlet,
    lon = long_outlet,
    state = state_outlet
  )




# 2. Aggregate daily data to annual data ---------------------------------------

## CONSTANTS ===================================================================
acceptable_missing_streamflow_days <- 10
minimum_entires_year <- 30 # at least 30 years of data required. Copied from HRS
min_run_length <- 10 # mainly for autocorrelation parameter estimation



## Yearly data =================================================================
yearly_precip <- daily_precip_raw |> # The entire precip data set is continuous
  pivot_longer(
    cols = !c(year:day),
    names_to = "gauge",
    values_to = "p_mm"
  ) |>
  summarise(
    p_mm = sum(p_mm),
    .by = c(year, gauge)
  ) |>
  arrange(gauge, year) 

yearly_streamflow <- daily_streamflow_raw |>
  pivot_longer(
    cols = !c(year:day),
    names_to = "gauge",
    values_to = "q_mm"
  ) |>
  mutate(is_na = is.na(q_mm)) |>
  summarise(
    q_mm = sum(q_mm, na.rm = TRUE), # add up everything. Ignore missing values
    q_na_count = sum(is_na),
    .by = c(year, gauge)
  ) |>
  mutate(
    # if there are more than acceptable_missing_streamflow_days missing values then set to NA
    q_mm = if_else(q_na_count >= acceptable_missing_streamflow_days, NA, q_mm) 
  ) |>
  arrange(gauge, year)


yearly_data <- yearly_precip |>
  right_join(yearly_streamflow, by = join_by(gauge, year)) 



## Continuous runs of data =====================================================
## Find the start and end of all continous streamflow measurements.
## Return a n x m matrix where m[1] = start_index, m[2] = end_index


gauge_continous_start_end <- function(single_gauge, data, min_run_length) {
  gauge_specific_data <- data |>
    filter(single_gauge == gauge)
  
  start_end_matrix <- continuous_run_start_end(gauge_specific_data$q_mm) # apply it to every gauge.
  
  start_end_matrix |>
    as_tibble() |>
    mutate(length = end_index - start_index + 1) |> # + 1 to make the input more easy to understand
    filter(length > min_run_length - 1) 
}


list_start_end_index <- map(
  .x = unique(yearly_data$gauge),
  .f = gauge_continous_start_end,
  data = yearly_data,
  min_run_length = min_run_length 
) 

names(list_start_end_index) <- paste0("gauge_", unique(yearly_data$gauge))

start_end_index <- do.call("rbind", list_start_end_index)

start_end_index <- start_end_index |>
  rownames_to_column(var = "gauge") |>
  as_tibble() |>
  mutate(gauge = str_remove(gauge, "\\.[1-9]")) |>
  mutate(gauge = str_remove(gauge, "gauge_")) # remove anything with "." followed by a number

write_csv(
  start_end_index,
  file = "./Data/Tidy/start_stop_index.csv"
)

gauges_in_start_end <- start_end_index |> pull(gauge) |> unique()

## Not all gauges have at a continuous period of at least min_run_length years
## Remove from yearly_data
## This is mainly to avoid having over-inflated autocorrelation values
yearly_data <- yearly_data |> 
  filter(gauge %in% gauges_in_start_end)

# 3. Find boxcox transform value for each gauge --------------------------------
gauge_info <- yearly_data |>
  summarise(
    bc_lambda = boxcox_lambda_generator(
      precipitation = p_mm, 
      streamflow = q_mm, 
      lambda_2 = 1 # lambda_2 + 1 to all streamflow to removes zeros
      ), 
    .by = gauge
  ) |>
  left_join( # Add chunks to gauge_info i.e., continuous runs
    catchment_information, 
    by = join_by(gauge)
  )




# 4. Convert streamflow (mm) into box-cox streamflow ---------------------------

## Function to convert streamflow (q_mm) into box-cox streamflow ===============
bc_q_generator <- function(gauge_id, a_priori_boxcox_lambda, yearly_data, lambda_2) {
  
  extracted_q_mm <- yearly_data |>
    filter(gauge == gauge_id) |>
    select(q_mm)
  
  if (any(extracted_q_mm[!is.na(extracted_q_mm)] <= 0) & (lambda_2 < 1)) {
    stop("Cannot transform value less than zero. Make lambda_2 >= 1")
  }
  
  boxcox_transform(extracted_q_mm, lambda = a_priori_boxcox_lambda, lambda_2 = lambda_2)
}



with_NA_bc_q <- map2(
  .x = gauge_info$gauge,
  .y = gauge_info$bc_lambda,
  .f = bc_q_generator,
  yearly_data = yearly_data,
  lambda_2 = 1
) |> 
  list_rbind()


names(with_NA_bc_q) <- "bc_q"


with_NA_yearly_data <- yearly_data |>
  add_column(
    with_NA_bc_q,
    .before = 5
  )


# 5. Add climate type to gauge_info using kgc ----------------------------------

## `LookupCZ` function provides the climate zone based on lon and lat 
## Relies on the climatezone dataframe
climatezones <- kgc::climatezones

## To use LookupCZ the data must be in |site_ID|lon|lat| format ================
formatted_gauge_information <- gauge_info |> 
  select(gauge, lat, lon) |> # mass overwrites dplyr select 
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


## Join to existing gauge_info =================================================
gauge_info <- gauge_info |> 
  left_join(
    climate_type_gauge_info,
    by = join_by(gauge)
  )


# 6. Save .csv  ----------------------------------------------------------------
write_csv(
  gauge_info, 
  paste0("./Data/Tidy/gauge_information_CAMELS.csv")
  )

write_csv(
  with_NA_yearly_data, 
  paste0("./Data/Tidy/with_NA_yearly_data_CAMELS.csv")
  )




