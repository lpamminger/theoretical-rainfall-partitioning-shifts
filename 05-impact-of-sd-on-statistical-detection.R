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
source("./Functions/synthetic_streamflow_model.R") 
source("./Functions/modified_stochastic_rainfall_generator.R")
source("./Functions/utility.R")



# Generate observations --------------------------------------------------------
## Generation parameters =======================================================
pre_shift_length_years <- 100 
post_shift_length_years <- 100
streamflow_parameter_multipliers <- c(1.1, 1.25, 1.5)

rainfall_parameters <- c("mean" = 1006,
                         "sd" = 221,
                         "auto" = 0.015,
                         "skew" = 0.19
                         )




# Statistical test function ----------------------------------------------------
p_value_generator_residuals <- function(parameter_position_index, streamflow_multiplier, post_rainfall_variation_years, normal_sd, detection_function){
  
  if (normal_sd) {
    streamflow_parameters <- c(
      "a0" = -4.1,
      "a1" = 0.017,
      "a2" = 0.16,
      "a3" = 2,
      "a4" = 0.014
    ) # Must be between -0.99 and 0.99
  } else if (!normal_sd) {
    streamflow_parameters <- c(
      "a0" = -4.1,
      "a1" = 0.017,
      "a2" = 0.16,
      "a3" = 0.5,
      "a4" = 0.014
    ) # Must be between -0.99 and 0.99
  }

  
  ## p_value_generator_residuals information ===================================
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
                                                         length_of_generated_rainfall = pre_shift_length_years,
                                                         random = TRUE)
  
  
  ### Change-rainfall ##########################################################
  change_rainfall <- modified_stochastic_rainfall_generator(parameter_vector = rainfall_parameters, 
                                                            length_of_generated_rainfall = post_shift_length_years,
                                                            random = TRUE)
  
  
  ### Apply synthetic_streamflow model #########################################

  
  ### Change-streamflow ########################################################
  ## Apply a single multiplier to a single streamflow parameter
  change_streamflow_parameters <- change_parameter_set_function(multiplier = streamflow_multiplier,
                                                                parameter_to_change = parameter_position_index, # the position of parameter in streamflow_parameters
                                                                control_parameter_set = streamflow_parameters)
  
  if (abs(change_streamflow_parameters[length(change_streamflow_parameters)]) >= 1) {stop("skewness cannot be greater than abs(1)")}
  

  
  change_synthetic_streamflow_model <- synthetic_streamflow_model(
                                        control_parameters = streamflow_parameters, 
                                        control_rainfall = pre_rainfall,
                                        random = TRUE
                                        )
                                      
  
  all_streamflow <- change_synthetic_streamflow_model(
                      change_parameters = change_streamflow_parameters,
                      change_rainfall = change_rainfall
                      )
  
  pre_streamflow <- all_streamflow[, 3]
  mean_only_pre_streamflow <- all_streamflow[, 2]
  change_streamflow <- all_streamflow[, 6]
  
  
  
  ### Residual calculation #####################################################
  # - fit a linear model to the mean_only_pre_streamflow and pre_rainfall
  # - use the linear model fitted to pre conditions to predict change streamflow using change rainfall
  # - change residauls = modelled_changed_streamflow - linear_model predicted streamflow using pre_conditions
  
  
  #### Find linear relationship between control rainfall and runoff ############
  linear_model_control <- lm(mean_only_pre_streamflow ~ pre_rainfall) # line
  
  pre_residuals <- pre_streamflow -  mean_only_pre_streamflow # pre-dot to pre-line distance
  
  pre_linear_model_coefficients <- c("intercept" = coef(linear_model_control)[1], 
                                     "slope" = coef(linear_model_control)[2])
  
  #### Apply pre_linear_model_coefficients to change_rainfall ##################
  #### to predict future streamflow 
  predicted_streamflow <- pre_linear_model_coefficients[1] + 
                            (pre_linear_model_coefficients[2] * change_rainfall)
  
  
  #### Find residual between change_streamflow and predicted_streamflow ########
  change_residuals <- change_streamflow - predicted_streamflow
  
  
  ### Apply test ###############################################################
  if (detection_function == "ks.test"){
    ks.test(pre_residuals, change_residuals[1:post_rainfall_variation_years])$p.value
    
  } else if (detection_function == "fligner.test"){
    # Fligner is special and wont work with two vectors
    # One vector contains the values and the other contains whether its pre/change
    modified_change_residuals <- change_residuals[1:post_rainfall_variation_years]
    
    combined_pre_change_residuals <- c(pre_residuals, 
                                        modified_change_residuals)
    
    groups <- c(rep("pre", time = length(pre_residuals)), 
                rep("change", time = length(modified_change_residuals))) 
    
    fligner.test(x = combined_pre_change_residuals, g = groups)$p.value
  }
}

# Find all combination of multiplier, parameters and years tested -------------- 
streamflow_parameters_index <- 3 #seq(from = 1, to = length(streamflow_parameters))
rainfall_variations <- seq(from = 10, to = post_shift_length_years)
low_or_high_sd <- c(TRUE, FALSE)

all_combinations_parameter_multi_rainfall <- expand.grid(
  streamflow_parameters_index,
  streamflow_parameter_multipliers,
  rainfall_variations,
  low_or_high_sd
)



# Apply model to all combinations and replicates -------------------------------
list_all_combinations_parameter_multi_rainfall <- list(all_combinations_parameter_multi_rainfall[,1], 
                                                       all_combinations_parameter_multi_rainfall[,2], 
                                                       all_combinations_parameter_multi_rainfall[,3],
                                                       all_combinations_parameter_multi_rainfall[,4]
                                                       )



## Define number of replicates =================================================
max_replicates <- 1000L

replicates <- seq(from = 1, to = max_replicates, by = 1)

repeat_p_value_generator_wrapper <- function(replicate, detection_function){
  
  p_values <- pmap_dbl(.l = list_all_combinations_parameter_multi_rainfall, 
                       .f = p_value_generator_residuals, 
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
                                       seed = 1L, 
                                       globals = TRUE),
                                     .progress = TRUE)

fligner_p_values_replicates <- future_map(.x = replicates, 
                                          .f = repeat_p_value_generator_wrapper,
                                          detection_function = "fligner.test",
                                          .options = furrr_options(
                                            seed = 1L, 
                                            globals = TRUE),
                                          .progress = TRUE)



# Summarise results in a dataframe ---------------------------------------------
making_summary_tibble <- function(surrogate_data_from_p_value_gen) {
  ## Surrogate data for joining
  names(surrogate_data_from_p_value_gen) <- paste("run", replicates, sep = "_")

  surrogate_p_value_tibble <- surrogate_data_from_p_value_gen |>
    as_tibble() |>
    mutate(surrogate_key = row_number())

  ## Main tibble
  summary_p_values <- as_tibble(all_combinations_parameter_multi_rainfall) |>
    rename(
      parameter = Var1,
      multiplier = Var2,
      years_post_change = Var3,
      low_or_high_sd = Var4
    ) |>
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


all_tests <- rbind(summary_ks_p_values, summary_fligner_p_values) |> 
  mutate(
    low_or_high_sd = if_else(low_or_high_sd, "Large", "Small")
    )



summary_all_tests <- all_tests |> 
  summarise(
    ave_p_value = mean(p_value),
    upper_p_value = quantile(p_value, 0.1), # can change
    lower_p_value = quantile(p_value, 0.9), # can change
    .by = c(ks_or_flig, parameter, multiplier, years_post_change, low_or_high_sd)
  )


# Plotting results -------------------------------------------------------------

## Combined graphs =============================================================
percentage_change <- (streamflow_parameter_multipliers * 100) - 100
multiplier_labs <- paste(paste0(percentage_change, "%"), "increase")
names(multiplier_labs) <- streamflow_parameter_multipliers

combined_low_high_sd_plot <- summary_all_tests |>
  ggplot(aes(x = years_post_change, y = ave_p_value, colour = low_or_high_sd)) +
  geom_line() +
  geom_ribbon(
    aes(
      x = years_post_change,
      ymin = lower_p_value,
      ymax = upper_p_value,
      fill = low_or_high_sd,
      colour = NULL
    ),
    alpha = 0.15
  ) +
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed") +
  labs(
    x = "Years of Post Change Data",
    y = "P-value",
    colour = "Magnitude of Control Standard Deviation",
    fill = "Magnitude of Control Standard Deviation"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(ks_or_flig ~ multiplier,
    labeller = labeller(
      multiplier = multiplier_labs
    )
  ) +
  theme(legend.position = "bottom")


single_low_high_sd_plot <- summary_all_tests |>
  filter(multiplier == 1.5) |>
  filter(ks_or_flig == "Kolmogorov-Smirnov") |>
  ggplot(aes(x = years_post_change, y = ave_p_value, colour = low_or_high_sd)) +
  geom_line() +
  geom_ribbon(
    aes(
      x = years_post_change,
      ymin = lower_p_value,
      ymax = upper_p_value,
      fill = low_or_high_sd,
      colour = NULL
    ),
    alpha = 0.15
  ) +
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed") +
  labs(
    x = "Years of Post Change Data",
    y = "P-value",
    colour = "Magnitude of Control Standard Deviation",
    fill = "Magnitude of Control Standard Deviation"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom")


# Save graphs ------------------------------------------------------------------

ggsave(paste0("./Graphs/how_low_high_sd_impacts_auto_detection", str_remove_all(Sys.Date(), "-"), ".pdf"),
       plot = combined_low_high_sd_plot,
       device = cairo_pdf,
       units = "mm",
       width = 210,
       height = 110)


ggsave(paste0("./Graphs/single_how_low_high_sd_impacts_auto_detection", str_remove_all(Sys.Date(), "-"), ".pdf"),
       plot = single_low_high_sd_plot,
       device = cairo_pdf,
       units = "mm",
       width = 210,
       height = 160)





