# detecting randomised change

# Clear environment and console ------------------------------------------------
cat("\014")

# Import libraries--------------------------------------------------------------
pacman::p_load(sn, moments, tidyverse, furrr, tictoc)

# Import functions -------------------------------------------------------------
source("./Functions/adjusting_parameters.R")
source("./Functions/utility.R")
source("./Functions/synthetic_streamflow_model.R") 
source("./Functions/modified_stochastic_rainfall_generator.R")
source("./Functions/utility.R")




pivot_streamflow_setup <- function(model_results) {
  ## get into |Output = bc_streamflow|Input = rainfall|Condition = control or change|
  ## so it has to be 3 cols and 200 rows

  model_results |> 
    rename(
      control_boxcoxstreamflow = control_boxcox_streamflow,
      change_boxcoxstreamflow = change_boxcox_streamflow
    ) |> 
    pivot_longer(
      cols = !replicate,
      names_to = c("control_or_change", ".value"),
      names_sep = "_" # there are two underscores
    ) |> 
    arrange(control_or_change)
}



generate_rainfall_runoff <- function(intercept_or_slope, replicate, intercept_multiplier, slope_multiplier, auto_rainfall_multiplier, pre_shift_length, post_shift_length, control_streamflow_parameters = NULL, control_rainfall_parameters = NULL) {
  
  # intercept_or_slope is a value of 1, 2, 3 or 4 
  
  ## Default change streamflow_parameters ======================================
  if (is.null(control_rainfall_parameters)) {
    control_rainfall_parameters <- c("mean" = 1006, 
                                     "sd" = 221, 
                                     "auto" = 0.015, 
                                     "skew" = 0.19
                                     )
  }
  
  ## Default control streamflow_parameters =====================================
  if (is.null(control_streamflow_parameters)) {
    control_streamflow_parameters <- c("a0" = -4.1, 
                                       "a1" = 0.017, 
                                       "a2" = 0.16, 
                                       "a3" = 2, 
                                       "a4" = 0.014
                                       )
  }
  
  
  ## Make change_streamflow_parameters =========================================
  change_streamflow_parameters <- change_parameter_set_function(
                                    multiplier = case_when(
                                      intercept_or_slope == 1 ~ intercept_multiplier,
                                      intercept_or_slope == 2 ~ slope_multiplier,
                                      .default = 1 # account for 3 here 
                                      ), 
                                    parameter_to_change = intercept_or_slope, # 3rd index is used here but it doesn't matter because the multipier = 1
                                    control_parameter_set = control_streamflow_parameters
                                    )
  
  change_rainfall_parameters <- change_parameter_set_function(
                                  multiplier = if_else(
                                    intercept_or_slope == 4,
                                    auto_rainfall_multiplier,
                                    1
                                  ),
                                  parameter_to_change = intercept_or_slope,
                                  control_parameter_set = control_rainfall_parameters
                                )
  
  ## Make rainfall =============================================================
  ### At the moment both control and change rainfall have the same parameters
  control_rainfall <- modified_stochastic_rainfall_generator(
                        parameter_vector = control_rainfall_parameters, 
                        length_of_generated_rainfall = pre_shift_length,
                        set_seed = FALSE
                        )
  

  change_rainfall <- modified_stochastic_rainfall_generator(
                       parameter_vector = change_rainfall_parameters,
                       length_of_generated_rainfall = post_shift_length,
                       set_seed = FALSE
                       )
  
  
  ## Make streamflow ===========================================================
  ### Control streamflow #######################################################
  streamflow_setup <- synthetic_streamflow_model(
                        control_parameters = control_streamflow_parameters,
                        control_rainfall = control_rainfall,
                        set_seed = FALSE
                        )
  
  
  ## Change streamflow ===========================================================
  model_results <- streamflow_setup(
                     change_parameters = change_streamflow_parameters,
                     change_rainfall = change_rainfall
                     ) 
  
  model_results <- model_results |> 
                     as_tibble() |> 
                     add_column(
                       replicate, 
                       .before = 1
                       )
  
  pivot_streamflow_setup(model_results)
}


# Running replicates -----------------------------------------------------------
## Run the replicate_intercept_slope_combinations for multiple intercept and slope multipliers


# Find the probablity of intercept or slope change for each replicate ----------
get_probability_of_intercept_change <- function(boxcox_streamflow, rainfall, control_or_change) {
  
  p_values <- lm(boxcox_streamflow ~ rainfall + control_or_change + (rainfall * control_or_change)) |> # or just intercept lm(boxcox_streamflow ~ rainfall + control_or_change) |> 
    summary() |> 
    coef()
  
  p_values[(nrow(p_values) - 1), ncol(p_values)] # or just intercept p_values[nrow(p_values), ncol(p_values)]
  
}


get_probability_of_slope_change <- function(boxcox_streamflow, rainfall, control_or_change) {
  
  p_values <- lm(boxcox_streamflow ~ rainfall + control_or_change + (rainfall * control_or_change)) |> 
    summary() |> 
    coef()
  
  p_values[nrow(p_values), ncol(p_values)]

}
  

change_multipliers_replicate_intercept_slope_combinations <- function(multiplier, intercept_slope_or_no_change, pre_shift_length, post_shift_length) {
  
  
  
  replicate_intercept_slope_combinations <- imap(
                                              .x = intercept_slope_or_no_change,
                                              .f = generate_rainfall_runoff,
                                              intercept_multiplier = multiplier, # keep same for plots
                                              slope_multiplier = multiplier, 
                                              auto_rainfall_multiplier = multiplier,
                                              pre_shift_length = pre_shift_length, 
                                              post_shift_length = post_shift_length
                                            ) |> 
                                            list_rbind()
  
  
  p_value_intercept_slope <- replicate_intercept_slope_combinations |> 
                               summarise(
                                 p_value_intercept_change = get_probability_of_intercept_change(
                                   boxcox_streamflow = boxcoxstreamflow,
                                   rainfall = rainfall,
                                   control_or_change = control_or_change
                                 ),
                                 p_value_slope_change = get_probability_of_slope_change(
                                   boxcox_streamflow = boxcoxstreamflow,
                                   rainfall = rainfall,
                                   control_or_change = control_or_change
                                 ),
                                 .by = replicate
                               )
  
  
  result <- p_value_intercept_slope |> 
              add_column(
                known_change = intercept_slope_or_no_change,
                .before = 1
              ) |> 
              mutate(
                known_change = case_when(
                  known_change == 1 ~ "intercept",
                  known_change == 2 ~ "slope",
                  known_change == 3 ~ "no_change",
                  known_change == 4 ~ "auto_rainfall"
                ),
                predicted_change = case_when(
                  (p_value_intercept_change < 0.05) & (p_value_slope_change < 0.05) ~ "both",
                  p_value_intercept_change < 0.05 ~ "intercept",
                  p_value_slope_change < 0.05 ~ "slope",
                  .default = "no_change"
                )
              )
  
    
  
  result |> 
    mutate(
      correct_identification = known_change == predicted_change
    ) |> 
    summarise(
      sum_known_change = n(),
      sum_correct_identification = sum(correct_identification),
      .by = known_change
    ) |> 
    mutate(
      percentage_correct = (sum_correct_identification / sum_known_change) * 100
    ) |>  
    select(known_change, percentage_correct) |> 
    add_column(
      multiplier = multiplier,
      .before = 1
    )
 
  
} 





## Number of replicates is dictated by the length of intercept_or_slope
REPLICATES <- 5000

intercept_slope_or_no_change <- sample(c(1, 2, 3), size = REPLICATES, replace = TRUE) # not testing 4

multipliers <- seq(from = 1.01, to = 2, by = 0.01)


plan(multisession, workers = availableCores())
identification_results_100_year <- future_map(
  .x = multipliers,
  .f = change_multipliers_replicate_intercept_slope_combinations,
  intercept_slope_or_no_change = intercept_slope_or_no_change,
  pre_shift_length = 100,
  post_shift_length = 100,
  .options = furrr_options(
    globals = TRUE,
    seed = TRUE
  ),
  .progress = TRUE
  ) |> 
  list_rbind() |> 
  add_column(
    shift_length = 100
  )


plan(multisession, workers = availableCores())
identification_results_20_year <- future_map(
  .x = multipliers,
  .f = change_multipliers_replicate_intercept_slope_combinations,
  intercept_slope_or_no_change = intercept_slope_or_no_change,
  pre_shift_length = 20,
  post_shift_length = 20,
  .options = furrr_options(
    globals = TRUE,
    seed = TRUE
  ),
  .progress = TRUE
  ) |> 
  list_rbind() |> 
  add_column(
    shift_length = 20
  )

identification_results <- rbind(identification_results_100_year, identification_results_20_year)



# Plot results -----------------------------------------------------------------
plot <- identification_results |> 
  filter(known_change != "auto_rainfall") |> # This doesn't do anything
  mutate(
    multiplier = (multiplier - 1) * 100,
    known_change = case_when(
      known_change == "intercept" ~ "Intercept",
      known_change == "slope" ~ "Slope",
      known_change == "no_change" ~ "No Change",
      known_change == "auto_rainfall" ~ "Rainfall Autocorrelation",
      .default = NA
    )
  ) |>
  ggplot(aes(
    x = multiplier,
    y = percentage_correct,
    colour = known_change,
    linetype = as.factor(shift_length)
  )) +
  geom_line(linewidth = 1) +
  labs(
    x = "Change in Parameter (%)",
    y = "Correct Change Identified (%)",
    colour = "Type of change",
    linetype = "Shift length (Years)"
  ) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )


ggsave(
  filename = paste0("correct_idenfitication_changes_", get_date(), ".pdf"),
  plot = plot,
  device = "pdf",
  path = "./Graphs",
  width = 210,
  height = 150,
  units = "mm"
)


