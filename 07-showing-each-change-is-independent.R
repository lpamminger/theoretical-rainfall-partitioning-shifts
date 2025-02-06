# Streamflow sensitivity analysis

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")
par(mfrow = c(1,1))


# Import libraries--------------------------------------------------------------
pacman::p_load(sn, moments, tidyverse, furrr, parallel)
# purrr and ggplot are part of the tidyverse
# stochastic_rainfall_generator function requires moments package
# ggplot2 and ggpubr for plotting (ggpubr is an extension of ggplot)
# synthetic_streamflow_model function requires sn package for skewed normal distribution



# Import functions -------------------------------------------------------------
source("./Functions/adjusting_parameters.R")
source("./Functions/modified_stochastic_rainfall_generator.R")
source("./Functions/synthetic_streamflow_model.R")
source("./Functions/utility.R")




## Making up control parameters ================================================
control_parameters <- c(
  "a0" = -4.1,
  "a1" = 0.017,
  "a2" = 0.16,
  "a3" = 2,
  "a4" = 0.014
  )

## Make up change parameters ===================================================
multipliers_for_control_parameters <- c(1.5, 1.3, 3, 2, 65) 

change_parameters <- imap(.x = multipliers_for_control_parameters, 
                          .f = change_parameter_set_function, 
                          control_parameter_set = control_parameters
                          )


repeat_rainfall_runoff <- function(change_streamflow_parameters, control_streamflow_parameters, skip) {
  # Inputs are change_parameters vector, control_parameters vector, length of generated rainfall, 
  
  
  # Generate control and change rainfall
  # Generate rainfall using a annual stochastic rainfall model -------------------
  rainfall <- modified_stochastic_rainfall_generator(
    parameter_vector = c(
      "mean" = 1006, 
      "sd" = 221, 
      "auto" = 0.015, 
      "skew" = 0.19), 
    length_of_generated_rainfall = 1000 + skip,
    set_seed = FALSE
    )

  
  # Generate control and change streamflow with UNIQUE rainfall
  change_synthetic_streamflow_model <- synthetic_streamflow_model(
    control_parameters = control_streamflow_parameters,
    control_rainfall = rainfall,
    set_seed = FALSE
    )
  
  
  rainfall_runoff_results <- change_synthetic_streamflow_model(
    change_parameters = change_streamflow_parameters,
    change_rainfall = rainfall 
    )
  
  return(rainfall_runoff_results)
}



# Repeat for each parameter
replicates_rainfall_runoff <- function(replicate_number, change_streamflow_parameters, control_streamflow_parameters, skip) {
  
  repeat_streamflow_parameter <- map(.x = change_streamflow_parameters,
                                     .f = repeat_rainfall_runoff,
                                     control_streamflow_parameters = control_streamflow_parameters,
                                     skip = skip
                                     )
  
  names_vector <- c("intercept",
                    "slope",
                    "autocorrelation",
                    "standard_deviation",
                    "skewness"
                    )
  
  
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

REPLICATES <- 5000 

replicate_results <- future_map(.x = seq(from = 1, to = REPLICATES, by = 1),
                                .f = replicates_rainfall_runoff,
                                change_streamflow_parameters = change_parameters,
                                control_streamflow_parameters = control_parameters,
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
## - Look at the statistical properties of the residuals around the line of best fit
###  * How does altering the standard deviation parameter impact the standard deviation around the residuals?
###  * how does altering the autocorrelation impact the skewness, autocorrelation and standard deviation of the residuals?



# Slope and intercept analysis -------------------------------------------------

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
    con_slo = mean(control_slope),
    .by = parameter_changed 
  ) |> 
  mutate(
    rel_diff_intercept = (res_intercept / con_int) * 100, 
    rel_diff_slope = (res_slope / con_slo) * 100
  )


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


# set constant binwidth for everything
hist_plot <- pivot_for_plotting |> 
               mutate(
                 intercept_or_slope = if_else(intercept_or_slope == "slope", "Fitted Slope", "Fitted Intercept"),
                 parameter_changed = case_when(
                   parameter_changed == "autocorrelation" ~ "Autocorrelation",
                   parameter_changed == "intercept" ~ "Intercept",
                   parameter_changed == "skewness" ~ "Skewness",
                   parameter_changed == "slope" ~ "Slope",
                   parameter_changed == "standard_deviation" ~ "Standard Deviation"
                   )
               ) |> 
               mutate(parameter_changed = factor( # factor used to set order of parameter_changed
                 parameter_changed, 
                 c("Intercept", 
                   "Slope", 
                   "Autocorrelation", 
                   "Standard Deviation", 
                   "Skewness"
                   )
                 )
               ) |> 
               ggplot(aes(x = intercept_or_slope_value)) + 
               geom_histogram(
                 fill = "grey", 
                 colour = "black", 
                 linewidth = 0.25,
                 position = "identity", 
                 bins = 50
                 ) +
               geom_vline(xintercept = 0, colour = "red") +
               labs(
                 x = "Residual (Change - Control)",
                 y = "Frequency"
                 ) +
               theme_bw() +
               theme(legend.title = element_blank()) +
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
  filename = paste0("./Graphs/intercept_slope_assessment_replicate_histogram_", REPLICATES, "_", "year_", nrow(summarise_intercept_slope)/5, "_", get_date(), ".pdf"),
  device = "pdf",
  width = 210,
  height = 190,
  units = "mm"
  )


