# Statistical detection

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")
par(mfrow = c(1,1))


# Import libraries--------------------------------------------------------------
pacman::p_load(sn, moments, tidyverse, furrr, parallel)


# Import functions -------------------------------------------------------------
source("./Functions/box_cox_transforms.R")
source("./Functions/adjusting_parameters.R")
source("./Functions/modified_stochastic_rainfall_generator.R")
source("./Functions/synthetic_streamflow_model.R")


# Generate observations --------------------------------------------------------
## Generation parameters =======================================================
pre_shift_length_years <- 100
post_shift_length_years <- 100
streamflow_parameter_multipliers <- c(1.1, 1.5, 2)

rainfall_parameters <- c(
  "mean" = 1006,
  "sd" = 221,
  "auto" = 0.015,
  "skew" = 0.19
)

streamflow_parameters <- c(
  "a0" = -4.1,
  "a1" = 0.017,
  "a2" = 0.16,
  "a3" = 2,
  "a4" = 0.014
) # Must be between -0.99 and 0.99


parameter_names <- names(streamflow_parameters)




# Statistical test function ----------------------------------------------------
p_value_generator <- function(parameter_position_index, streamflow_multiplier, post_rainfall_variation_years, detection_function){
  
  
  ## p_value_generator information =====================================
  # This functions relies on global variables:
  ## - the rainfall/streamflow generators
  ## - pre/post shift lengths
  ## - the streamflow and rainfall parameter vectors
  
  # Inputs:
  ## - parameter_position_index - c(a0, a1, a2, a3, a4) = c(1, 2, 3, 4, 5)
  ## - streamflow_multiplier - multiplier to parameter
  ## - post_rainfall_variation_year - years
  
  
  ### Pre-rainfall #############################################################
  pre_rainfall <- modified_stochastic_rainfall_generator(parameter_vector = rainfall_parameters, 
                                                         length_of_generated_rainfall = pre_shift_length_years
                                                         )
  
  
  ### Change-rainfall ##########################################################
  change_rainfall <- modified_stochastic_rainfall_generator(parameter_vector = rainfall_parameters, 
                                                            length_of_generated_rainfall = post_shift_length_years
                                                            )
  
  
  ### Apply synthetic streamflow model #########################################
  
  ### Change-streamflow ########################################################
  ## Apply a single multiplier to a single streamflow parameter
  change_streamflow_parameters <- change_parameter_set_function(multiplier = streamflow_multiplier,
                                                                parameter_to_change = parameter_position_index, # the position of parameter in streamflow_parameters
                                                                control_parameter_set = streamflow_parameters
                                                                )
  
  
  change_synthetic_streamflow_model <- synthetic_streamflow_model(
                                         control_parameters = streamflow_parameters, 
                                         control_rainfall = pre_rainfall,
                                         random = NULL
                                         )
  
  
  all_streamflow <- change_synthetic_streamflow_model(
                      change_parameters = change_streamflow_parameters,
                      change_rainfall = change_rainfall
                      )
  
  pre_streamflow <- all_streamflow[, 2]
  
  change_streamflow <- all_streamflow[, 4]
  
  
  
  ### Apply test ###############################################################
  if (detection_function == "ks.test"){
    ks.test(pre_streamflow, change_streamflow[1:post_rainfall_variation_years])$p.value # The first 
    
  } else if (detection_function == "fligner.test"){
    # Fligner is special and wont work with two vectors
    # One vector contains the values and the other contains whether its pre/change
    modified_change_streamflow <- change_streamflow[1:post_rainfall_variation_years]
    
    combined_pre_change_streamflow <- c(pre_streamflow, 
                                        modified_change_streamflow)
    
    groups <- c(rep("pre", time = length(pre_streamflow)), 
                rep("change", time = length(modified_change_streamflow))) 
    
    fligner.test(x = combined_pre_change_streamflow, g = groups)$p.value
  } else {
    stop("Function name not found")
  }
}


# For paper special case for autocorrelation and standard dev ------------------
# What I want to do:
# - apply trend test to autocorrelation term
# - change the standard deviation
# - see how changing the standard deviation impact the likelihood of detection of autocorrelation
# sd is fixed default

## Objective is to see how changing autocorrelation impacts sd
wrapper_special_case <- function(post_rainfall_variation_years, multiplier, detection_function) {
  
  p_value_generator(parameter_position_index = 3, # 3 = autocorrelation
                    streamflow_multiplier = multiplier, # WRONG we want to change sd
                    post_rainfall_variation_years = post_rainfall_variation_years,
                    detection_function = detection_function
                    )
}


pre_allocate_replicate <- matrix(numeric(length = 91 * 100), nrow = 91, ncol = 100)

for (i in 1:100) {
  single_replicate <- map_dbl(.x = seq(from = 10, to = post_shift_length_years),
                              .f = wrapper_special_case,
                              multiplier = 2,
                              detection_function = "ks.test"
                              )
  
  pre_allocate_replicate[,i] <- single_replicate
}  






# Find all combination of multiplier, parameters and years tested -------------- 
streamflow_parameters_index <- seq(from = 1, to = length(streamflow_parameters))
rainfall_variations <- seq(from = 10, to = post_shift_length_years)

all_combinations_parameter_multi_rainfall <- as.matrix(
                                               expand.grid(
                                                 streamflow_parameters_index, 
                                                 streamflow_parameter_multipliers, 
                                                 rainfall_variations
                                                 )
                                               )


# Apply model to all combinations and replicates -------------------------------
list_all_combinations_parameter_multi_rainfall <- list(all_combinations_parameter_multi_rainfall[,1], 
                                                       all_combinations_parameter_multi_rainfall[,2], 
                                                       all_combinations_parameter_multi_rainfall[,3]
                                                       )

## Define number of replicates =================================================
max_replicates <- 1000L

replicates <- seq(from = 1, to = max_replicates, by = 1)

repeat_p_value_generator_wrapper <- function(replicate, detection_function){
  
  p_values <- pmap_dbl(.l = list_all_combinations_parameter_multi_rainfall, 
                          .f = p_value_generator, 
                          detection_function = detection_function,
                          .progress = FALSE) 
  
}

## Set up parallel processing using future_map =================================
plan(multisession, workers = availableCores())


## Run replicates in parallel ==================================================
ks_p_values_replicates <- future_map(.x = replicates,
                                     .f = repeat_p_value_generator_wrapper,
                                     detection_function = "ks.test",
                                     .options = furrr_options(
                                       seed = NULL, # Unless seed = FALSE/NULL it will generate the exact same sequence of random numbers 
                                       globals = TRUE),
                                     .progress = TRUE)

plan(multisession, workers = availableCores())


fligner_p_values_replicates <- future_map(.x = replicates,
                                          .f = repeat_p_value_generator_wrapper,
                                          detection_function = "fligner.test",
                                          .options = furrr_options(
                                            seed = NULL, # Unless seed = FALSE/NULL it will generate the exact same sequence of random numbers 
                                            globals = TRUE),
                                          .progress = TRUE)




# Summarise results in a dataframe ---------------------------------------------
making_summary_tibble <- function(surrogate_data_from_p_value_gen){
  
  ## Surrogate data for joining
  names(surrogate_data_from_p_value_gen) <- paste("run", replicates, sep = "_")
  
  surrogate_p_value_tibble <- surrogate_data_from_p_value_gen |> 
                                as_tibble() |> 
                                mutate(surrogate_key = row_number())
  
  ## Main tibble
  summary_p_values <- as_tibble(all_combinations_parameter_multi_rainfall) |> 
                        rename(parameter = Var1,
                               multiplier = Var2,
                               years_post_change = Var3) |> 
                        mutate(
                          parameter = case_when(
                            parameter == 1 ~ "a0", 
                            parameter == 2 ~ "a1",
                            parameter == 3 ~ "a2",
                            parameter == 4 ~ "a3",
                            parameter == 5 ~ "a4",
                          )
                        ) |> 
                        mutate(surrogate_key = row_number()) |> 
                        left_join(surrogate_p_value_tibble, join_by(surrogate_key)) |> 
                        select(!surrogate_key) |> 
                        pivot_longer(
                          cols = starts_with("run"),
                          names_to = "run",
                          values_to = "p_value"
                        ) 
}


# Summary p_values =============================================================
summary_ks_p_values <- making_summary_tibble(ks_p_values_replicates) |>
  add_column("ks_or_flig" = "Kolmogorov-Smirnov", .before = 1)

summary_fligner_p_values <- making_summary_tibble(fligner_p_values_replicates) |>
  add_column("ks_or_flig" = "Fligner-Killeen", .before = 1)


all_tests <- rbind(summary_ks_p_values, summary_fligner_p_values)

summary_all_tests <- all_tests |>
  summarise(
    ave_p_value = mean(p_value),
    upper_p_value = quantile(p_value, 0.1), # can change
    lower_p_value = quantile(p_value, 0.9), # can change
    .by = c(ks_or_flig, parameter, multiplier, years_post_change)
  )
# Plotting results -------------------------------------------------------------
parameter_labs <- c("Intercept", "Slope", "Autocorrelation", "Standard Deviation", "Skewness") # paste("Parameter", parameter_names)
names(parameter_labs) <- parameter_names

percentage_change <- (streamflow_parameter_multipliers * 100) - 100
multiplier_labs <- paste(paste0(percentage_change, "%"), "increase")
names(multiplier_labs) <- streamflow_parameter_multipliers


## Combined ks and fligner
combined_residual_detection_plot <- summary_all_tests |>
  ggplot(aes(x = years_post_change, y = ave_p_value, colour = ks_or_flig)) +
  geom_line() +
  geom_ribbon(
    aes(
      x = years_post_change,
      ymin = lower_p_value,
      ymax = upper_p_value,
      fill = ks_or_flig,
      colour = NULL
    ),
    alpha = 0.15
  ) +
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed") +
  labs(
    x = "Years of Post Change Data",
    y = "P-value",
    colour = "Statistical Test Used",
    fill = "Statistical Test Used"
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(parameter ~ multiplier,
    labeller = labeller(
      parameter = parameter_labs,
      multiplier = multiplier_labs
    )
  ) +
  theme(legend.position = "bottom")

ggsave(paste0("./Graphs/combined_streamflow_detection_", str_remove_all(Sys.Date(), "-"), ".pdf"),
  plot = combined_residual_detection_plot,
  device = cairo_pdf,
  units = "mm",
  width = 190,
  height = 230
)



