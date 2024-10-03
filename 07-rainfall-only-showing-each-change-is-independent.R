# Streamflow sensitivity analysis

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")


# Import libraries--------------------------------------------------------------
pacman::p_load(sn, moments, tidyverse, furrr, parallel)
# purrr and ggplot are part of the tidyverse
# stochastic_rainfall_generator function requires moments package
# ggplot2 and ggpubr for plotting (ggpubr is an extension of ggplot)
# synthetic_streamflow_model function requires sn package for skewed normal distribution



# Import functions -------------------------------------------------------------
source("./Functions/adjusting_parameters.R")



## Making up control parameters ================================================
control_streamflow_parameters <- c(
  "a0" = -4.1,
  "a1" = 0.017,
  "a2" = 0.16,
  "a3" = 2,
  "a4" = 0.014
)


control_rainfall_parameters <- c("mean" = 1006, 
                                 "sd" = 221, 
                                 "auto" = 0.015, 
                                 "skew" = 0.19
                                 )

## Make up change parameters ===================================================
multipliers_for_control_rainfall_parameters <- c(0.8, 1.3, 5, 5) 

change_rainfall_parameters <- imap(.x = multipliers_for_control_rainfall_parameters, 
                                   .f = change_parameter_set_function, 
                                   control_parameter_set = control_rainfall_parameters
                                   )



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










repeat_rainfall_runoff <- function(change_rainfall_parameters, control_rainfall_parameters, control_streamflow_parameters, skip) {
  # Inputs are change_parameters vector, control_parameters vector, length of generated rainfall, 
  
  
  # Generate control and change rainfall
  # Generate rainfall using a annual stochastic rainfall model -------------------
  control_rainfall <- modified_stochastic_rainfall_generator(
    parameter_vector = control_rainfall_parameters, 
    length_of_generated_rainfall = 1000 + skip
    )
  
  change_rainfall <- modified_stochastic_rainfall_generator(
    parameter_vector = change_rainfall_parameters, 
    length_of_generated_rainfall = 1000 + skip
  )

  
  # Generate control and change streamflow with UNIQUE rainfall
  change_synthetic_streamflow_model <- synthetic_streamflow_model(
    control_parameters = control_streamflow_parameters,
    control_rainfall = control_rainfall
    )
  
  
  rainfall_runoff_results <- change_synthetic_streamflow_model(
    change_parameters = control_streamflow_parameters,
    change_rainfall = change_rainfall 
    )
  
  return(rainfall_runoff_results)
}



# Repeat for each parameter
replicates_rainfall_runoff <- function(replicate_number, change_rainfall_parameters, control_rainfall_parameters, control_streamflow_parameters, skip) {
  
  repeat_streamflow_parameter <- map(.x = change_rainfall_parameters,
                                     .f = repeat_rainfall_runoff,
                                     control_rainfall_parameters = control_rainfall_parameters,
                                     control_streamflow_parameters = control_streamflow_parameters,
                                     skip = skip)
  
  names_vector <- names(control_rainfall_parameters)
  
  
  for (item in seq_along(repeat_streamflow_parameter)){
    repeat_streamflow_parameter[[item]] <- cbind(
                                             as_tibble(repeat_streamflow_parameter[[item]][(skip + 1):nrow(repeat_streamflow_parameter[[item]]),]), 
                                             names_vector[item]
                                             )
  }
  
  summary_tibble <- do.call(rbind, repeat_streamflow_parameter)
  
  summary_tibble <- summary_tibble |> 
    rename(
      parameter_changed = `names_vector[item]`
    ) |> 
    add_column(replicate = {{  replicate_number  }}, .before = 1) 
 
  return(summary_tibble)
   
}




## Run replicates in parallel ==================================================

# Replicates and random = TRUE should only be used for sd and skewness?


plan(multisession, workers = availableCores())

REPLICATES <- 5000 # 5000

replicate_results <- future_map(.x = seq(from = 1, to = REPLICATES, by = 1),
                                .f = replicates_rainfall_runoff,
                                change_rainfall_parameters = change_rainfall_parameters,
                                control_rainfall_parameters = control_rainfall_parameters,
                                control_streamflow_parameters = control_streamflow_parameters,
                                skip = 2,
                                .options = furrr_options(
                                  seed = NULL, # I don't want it to generate the exact same sequence of random numbers every time
                                  globals = TRUE),
                                .progress = TRUE
                                )

replicate_results <- replicate_results |> list_rbind()



# Analysis ---------------------------------------------------------------------
# Objectives:
## - Show the control and change slope and intercept trend to zero with the number of replicates


# Compare mean rainfall for each of the changes --------------------------------
mean_rainfall_check <- replicate_results |> 
                        summarise(
                          mean_control_rainfall = mean(control_rainfall),
                          mean_change_rainfall = mean(change_rainfall),
                          .by = parameter_changed
                        )






# rainfall-runoff slope analysis -----------------------------------------------
get_intercept <- function(rainfall, streamflow) {
  coef(lm(streamflow ~ rainfall))[1]
}

get_slope <- function(rainfall, streamflow) {
  coef(lm(streamflow ~ rainfall))[2]
}


summarise_intercept_slope <- replicate_results |> 
                               summarise(
                                 control_intercept = get_intercept(control_rainfall, control_boxcox_streamflow),
                                 change_intercept = get_intercept(change_rainfall, change_boxcox_streamflow),
                                 control_slope = get_slope(control_rainfall, control_boxcox_streamflow),
                                 change_slope = get_slope(change_rainfall, change_boxcox_streamflow),
                                 .by = c(parameter_changed, replicate)
                                 ) |> 
                               mutate(
                                 intercept_residual = change_intercept - control_intercept,
                                 slope_residual = change_slope - control_slope
                                 ) |> 
                               arrange(parameter_changed) 

some_stats <- summarise_intercept_slope |> 
                summarise(
                  res_intercept = mean(intercept_residual),
                  res_slope = mean(slope_residual),
                  con_int = mean(control_intercept),
                  cha_int = mean(change_intercept),
                  con_slo = mean(control_slope),
                  cha_slo = mean(change_slope),
                  .by = parameter_changed 
                  ) |> 
                mutate(
                  rel_diff_intercept = (res_intercept / con_int) * 100, 
                  rel_diff_slope = (res_slope / con_slo) * 100
                )


# t-test on residuals ---------------------------------------------------------
t_test_results <- summarise_intercept_slope |> 
  summarise(
    t_test_intercept = t.test(
      intercept_residual, 
      alternative = "two.sided",
      mu = 0
    )$p.value,
    t_test_slope = t.test(
      slope_residual, 
      alternative = "two.sided",
      mu = 0
    )$p.value,
    .by = parameter_changed
  )

# relative change = final - initial / initial

pivot_for_plotting <- summarise_intercept_slope |> 
                        select(!starts_with("control")) |> 
                        select(!starts_with("change")) |> 
                        pivot_longer(
                          cols = c(intercept_residual, slope_residual),
                          names_to = "intercept_or_slope",
                          values_to = "intercept_or_slope_value" 
                          ) |> 
                        mutate(
                          intercept_or_slope = str_remove_all(intercept_or_slope, "_residual")
                          )
  

hist_plot <- pivot_for_plotting |> 
              mutate(
                intercept_or_slope = if_else(intercept_or_slope == "slope", "Fitted Slope", "Fitted Intercept"),
                parameter_changed = case_when(
                  parameter_changed == "mean" ~ "Mean",
                  parameter_changed == "auto" ~ "Autocorrelation",
                  parameter_changed == "skew" ~ "Skewness",
                  parameter_changed == "sd" ~ "Standard Deviation",
                )
              ) |> 
              mutate(parameter_changed = factor( # factor used to set order of parameter_changed
                parameter_changed, 
                  c("Mean", 
                    "Standard Deviation", 
                    "Autocorrelation", 
                    "Skewness"
                    )
                  )
                ) |>
              ggplot(aes(x = intercept_or_slope_value)) + 
              geom_histogram(
                fill = "grey", 
                colour = "black", 
                linewidth = 0.2,
                position = "identity", 
                bins = 50
                ) +
              geom_vline(xintercept = 0, colour = "red") +
              labs(
                x = "Residual (Change - Control)",
                y = "Frequency"
                ) +
              theme_bw() +
              facet_grid(
                parameter_changed ~ intercept_or_slope,
                scales = "free"
                ) +
              theme(
                strip.text = element_text(size = 10), 
                panel.grid.minor = element_blank(),
                axis.title = element_text(size = 12)
              )


ggsave(
  plot = hist_plot, 
  filename = paste0("./Graphs/rainfall_only_intercept_slope_assessment_replicate_histogram_", REPLICATES, "_", "year_", nrow(summarise_intercept_slope)/4, "_", str_remove_all(Sys.Date(), "-"), ".pdf"),
  device = "pdf",
  width = 210,
  height = 160,
  units = "mm"
  )


