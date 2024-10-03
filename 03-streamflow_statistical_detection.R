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


# Generate observations --------------------------------------------------------
## Generation parameters =======================================================
pre_shift_length_years <- 100 
post_shift_length_years <- 100
streamflow_parameter_multipliers <- c(1.1, 1.5, 2)

rainfall_parameters <- c("mean" = 1006,
                         "sd" = 221,
                         "auto" = 0.015,
                         "skew" = 0.19)

streamflow_parameters <- c("a0" = -4.1,
                           "a1" = 0.017,
                           "a2" = 0.16,
                           "a3" = 2,
                           "a4" = 0.014) # Must be between -0.99 and 0.99


parameter_names <- names(streamflow_parameters)



# Altered modified stochastic rainfall generate - remove the set.seed ----------
modified_stochastic_rainfall_generator <- function(parameter_vector, length_of_generated_rainfall){
  
  # ============================================================================
  # Obtain statistical properties of observed rainfall
  #browser() # for debugging
  
  ## Get mean (x_bar) and std (s) of obs. annual rainfall
  x_bar <- parameter_vector[1]
  s <- parameter_vector[2]
  
  ## Get lag-1 autocorrelation coefficent (r) from obs. rainfall
  r <- parameter_vector[3]
  
  ## Get the skewness (gamma) from obs. rainfall
  gamma <- parameter_vector[4]
  
  ## Get coefficient of skewness (gamma_e)
  gamma_e <- ((1 - r^3)  / ((1 - r^2)^(3/2))) * gamma
  
  # ============================================================================
  
  
  # ============================================================================
  # Generating rainfall
  
  eta <- rnorm(length_of_generated_rainfall, mean = 0, sd = 1)
  
  ## Generate random number using eta_t and gamma_e (epsilon)
  epsilon <- (2 / gamma_e) * (((1 + ((gamma_e * eta) / 6) - (gamma_e^2 / 36))^3) - 1)
  
  ## Generate standardised annual rainfall (X)
  X <- numeric(length_of_generated_rainfall)
  X_prev <- 0 # assume X_prev = zero
  
  for (t in 1:length_of_generated_rainfall){
    
    X[t] <- (r * X_prev) + (sqrt(1 - (r^2)) * epsilon[t])
    X_prev <- X[t]
    
  }
  
  
  ## Annual rainfall amount (generated_rainfall)
  generated_rainfall <- x_bar + (s * X)
  # ============================================================================
  
  return(generated_rainfall)
}



# Altered synthetic streamflow model - remove set.seed -------------------------
synthetic_streamflow_model <- function(control_parameters, control_rainfall) {
  
  force(control_parameters)
  force(control_rainfall)
  
  # get control_parameters
  control_intercept <- control_parameters[1]
  control_slope <- control_parameters[2]
  control_autocorr <- control_parameters[3]
  control_stand_dev <- control_parameters[4]
  control_skew <- control_parameters[5]
  
  # control decay value
  decay_value <- control_intercept + (control_slope * mean(control_rainfall))
  
  # set previous value
  control_boxcox_streamflow_previous <- decay_value
  
  control_boxcox_streamflow <- numeric(length = length(control_rainfall))
  
  # Do I care about the mean-boxcox streamflow? not including it may change some things
  
  # get standard deviation and skewness to the skewed normal parameters
  control_skewed_normal_parameters <- to_skewed_normal_parameters(
    user_mean = 0,
    user_std = control_stand_dev,
    user_skewness = control_skew)
  
  # Run the model
  for (time_step in seq_along(control_rainfall)) {
    
    control_boxcox_streamflow[time_step] <- control_intercept + 
      (control_slope * control_rainfall[time_step]) + 
      (control_autocorr * (control_boxcox_streamflow_previous - decay_value)) +
      rsn(n = 1, xi = control_skewed_normal_parameters[1], omega = control_skewed_normal_parameters[2], alpha = control_skewed_normal_parameters[3])
    
    control_boxcox_streamflow_previous <- control_boxcox_streamflow[time_step]
    
  }
  
  
  
  function(change_parameters, change_rainfall) {
    
    #browser() # for debugging
    
    force(change_parameters)
    force(change_rainfall)
    
    
    # get change_parameters
    change_intercept <- change_parameters[1]
    change_slope <- change_parameters[2]
    change_autocorr <- change_parameters[3]
    change_stand_dev <- change_parameters[4]
    change_skew <- change_parameters[5]
    
    # change decay value
    decay_value <- change_intercept + (change_slope * mean(change_rainfall))
    
    # set previous value
    change_boxcox_streamflow_previous <- decay_value
    
    change_boxcox_streamflow <- numeric(length = length(change_rainfall))
    
    # Do I care about the mean-boxcox streamflow? not including it may change some things
    
    # get standard deviation and skewness to the skewed normal parameters
    change_skewed_normal_parameters <- to_skewed_normal_parameters(user_mean = 0,
                                                                   user_std = change_stand_dev,
                                                                   user_skewness = change_skew)
    
    # Run the model
    for (time_step in seq_along(change_rainfall)) {
      
      change_boxcox_streamflow[time_step] <- change_intercept + 
        (change_slope * change_rainfall[time_step]) + 
        (change_autocorr * (change_boxcox_streamflow_previous - decay_value)) +
        rsn(n = 1, xi = change_skewed_normal_parameters[1], omega = change_skewed_normal_parameters[2], alpha = change_skewed_normal_parameters[3])
      
      change_boxcox_streamflow_previous <- change_boxcox_streamflow[time_step]
      
    }
    
    
    # I want the function to return the mean and actual boxcox streamflow (as a tibble) |control_mean|control_actual|change_mean|change_actual
    result <- cbind(control_rainfall, control_boxcox_streamflow, change_rainfall, change_boxcox_streamflow)
    return(result)
    
  }
}





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
                                                                control_parameter_set = streamflow_parameters)
  
  
  change_synthetic_streamflow_model <- synthetic_streamflow_model(
                                         control_parameters = streamflow_parameters, 
                                         control_rainfall = pre_rainfall
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



# Plotting results -------------------------------------------------------------

## Plotting function ===========================================================
ggplotting_p_value_facet <- function(summary_p_values){
  
  ### Manipulate data to get the mean, upper and lower percentiles #############
  plotting_summary_p_values <- summary_p_values |> 
                                    summarise(
                                      ave_p_value = mean(p_value),
                                      upper_p_value = quantile(p_value, 0.1), # can change
                                      lower_p_value = quantile(p_value, 0.9), # can change
                                      .by = c(parameter, multiplier, years_post_change)
                                    )
  
  
  ### Labelling facet bits #####################################################
  parameter_labs <- c("Intercept", "Slope", "Autocorrelation", "Standard Deviation", "Skewness") #paste("Parameter", parameter_names)
  names(parameter_labs) <- parameter_names
  
  percentage_change <- (streamflow_parameter_multipliers * 100) - 100
  multiplier_labs <- paste(paste0(percentage_change, "%"), "increase")
  names(multiplier_labs) <- streamflow_parameter_multipliers
  
  
  ### Actual plot ##############################################################
  ggplot_results <- plotting_summary_p_values |> 
    ggplot(aes(x = years_post_change, y = ave_p_value)) +
    geom_line() + 
    geom_ribbon(
      aes(x = years_post_change, 
          ymin = lower_p_value,
          ymax = upper_p_value),
      alpha = 0.25) +
    geom_hline(yintercept = 0.05, 
               colour = "red",
               linetype = "dashed") +
    labs(x = "Years of data post change",
         y = "P-value significance") +
    theme_bw() +
    facet_grid(parameter ~ multiplier,
               labeller = labeller(
                 parameter = parameter_labs,
                 multiplier = multiplier_labs
               ))
  
}


## Results =====================================================================
ggplot_results_ks <- ggplotting_p_value_facet(summary_ks_p_values)
ggplot_results_fligner <- ggplotting_p_value_facet(summary_fligner_p_values)


# Saving -----------------------------------------------------------------------
ggsave(paste0("./Graphs/detection_test_ks_", str_remove_all(Sys.Date(), "-"), ".pdf"),
       plot = ggplot_results_ks,
       device = cairo_pdf,
       units = "mm",
       width = 190,
       height = 230)


ggsave(paste0("./Graphs/detection_test_fligner_", str_remove_all(Sys.Date(), "-"), ".pdf"),
       plot = ggplot_results_fligner,
       device = cairo_pdf,
       units = "mm",
       width = 190,
       height = 230)



