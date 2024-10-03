# Rainfall sensitivity analysis

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")
par(mfrow = c(1,1))


# Import libraries--------------------------------------------------------------
pacman::p_load(sn, moments, tidyverse)


# Import functions -------------------------------------------------------------
source("./Functions/box_cox_transforms.R")
source("./Functions/adjusting_parameters.R")
source("./Functions/synthetic_streamflow_model.R") 
source("./Functions/modified_stochastic_rainfall_generator.R")
source("./Functions/utility.R")


# Length of timeseries ---------------------------------------------------------
skip <- 2
user_time_series_length <- 100 # the length of rainfall determines the length of streamflow
time_series_length <- user_time_series_length + skip # we remove the first two values as the model spins up


## Control parameters are based on sample rainfall =============================
control_parameters <- c("mean" = 1006, 
                        "sd" = 221, 
                        "auto" = 0.015, 
                        "skew" = 0.19
                        )



## Setting up parameter sets list ==============================================
multipliers_for_control_parameters <- c(0.8, 1.3, 5, 5)

change_parameters <- imap(.x = multipliers_for_control_parameters, 
                          .f = change_parameter_set_function, 
                          control_parameter_set = control_parameters
                          )

parameter_list <- c(list(control_parameters), change_parameters)

names(parameter_list) <- c("control",
                           "mean",
                           "standard_deviation",
                           "autocorrelation",
                           "skewness"
                           )


## Put the parameters through the stochastic rainfall generator ============
rainfall_list <- map(.x = parameter_list,
                     .f = modified_stochastic_rainfall_generator,
                     length_of_generated_rainfall = time_series_length,
                     random = FALSE
                     )



# Generate streamflow using the generated rainfall -----------------------------
streamflow_parameters <- c(-4.1, 0.017, 0.16, 2, 0.014)


## Apply model =================================================================

### Part 1 of function factory - store the control values
change_synthetic_streamflow_model <- synthetic_streamflow_model(
                                       control_parameters = streamflow_parameters, 
                                       control_rainfall = rainfall_list$control,
                                       random = FALSE
                                       )



## Make the precipitation the first argument so I can loop over ================
modified_change_synthetic_streamflow_model <- function(change_rainfall, change_parameters) {
  change_synthetic_streamflow_model(change_parameters = change_parameters,
                                    change_rainfall = change_rainfall)
}


### Part 2 of function factory - add the change values
boxcox_streamflow_list <- map(.x = rainfall_list[-1], # exclude control
                              .f = modified_change_synthetic_streamflow_model, 
                              change_parameters = streamflow_parameters
                              )


# add another column
for (item in seq_along(boxcox_streamflow_list)){
  boxcox_streamflow_list[[item]] <- cbind(
                                      as_tibble(boxcox_streamflow_list[[item]][(skip + 1):nrow(boxcox_streamflow_list[[item]]),]), 
                                      names(parameter_list[-1])[item], 
                                      seq(from = 1, to = nrow(boxcox_streamflow_list[[item]][(skip + 1):nrow(boxcox_streamflow_list[[item]]),]))
                                      )
}

boxcox_streamflow <- do.call("rbind", boxcox_streamflow_list)

boxcox_streamflow <- boxcox_streamflow |> 
                       rename(
                         parameter = `names(parameter_list[-1])[item]`,
                         time = `seq(from = 1, to = nrow(boxcox_streamflow_list[[item]][(skip + `
                       ) |> 
                       remove_rownames() |> 
                       relocate(
                         c(parameter, time),
                         .before = 1
                         )



# Get in tidy form
## I do not know who to pivot longer into this format --> parameter|control_or_change|rainfall|boxcox_mean_streamflow|boxcox_streamflow|streamflow
## So I will use multiple steps

get_boxcox_streamflow_into_tidy <- function(name_1, name_2, value_to_name, data) {
  
  data |> 
    select(parameter, time, {{ name_1 }}, {{ name_2 }}) |> 
    pivot_longer(
      cols = !c(parameter, time),
      names_to = "control_or_change",
      values_to = value_to_name
    ) |> 
    mutate(
      control_or_change = str_extract(control_or_change, "control|change")
    )
}

name_1 <- boxcox_streamflow |> 
  select(starts_with("control")) |> 
  names()

name_2 <- boxcox_streamflow |> 
  select(starts_with("change")) |> 
  names()


value_to_name <- c("rainfall", "boxcox_streamflow")

tidy_boxcox_streamflow <- pmap(.l = list(name_1, name_2, value_to_name),
                               .f = get_boxcox_streamflow_into_tidy,
                               data = boxcox_streamflow)

tidy_boxcox_streamflow <- purrr::reduce(tidy_boxcox_streamflow, left_join, by = join_by(parameter, time, control_or_change))


### Convert from boxcox to real space ##########################################
#### When talking the average of hState and CAMELS catchments is around 0.3
boxcox_lambda <- 0.3 #o.g. = 0.3

tidy_boxcox_streamflow <- tidy_boxcox_streamflow |> 
                            mutate(
                              streamflow = BCTransformInverse(boxcox_streamflow, lambda = boxcox_lambda)
                            )                             
                            



# Statistical tests ------------------------------------------------------------
## Stats table =================================================================

# Are the differences in these values due to the random number set? Test using replicates...

statistical_properties <- tidy_boxcox_streamflow |> 
  summarise(
    intercept = coef(lm(boxcox_streamflow ~ rainfall))[1],
    slope = coef(lm(boxcox_streamflow ~ rainfall))[2],
    .by = c(parameter, control_or_change)
  ) |> 
  arrange()



# Plotting ---------------------------------------------------------------------

# force the order of items
tidy_boxcox_streamflow$parameter <- factor(
  tidy_boxcox_streamflow$parameter,
  levels = c(
    "mean",
    "standard_deviation",
    "autocorrelation",
    "skewness"
  )
)

tidy_boxcox_streamflow <- tidy_boxcox_streamflow |>
  mutate(
    control_or_change = if_else(control_or_change == "control", "Control", "Change")
  )


## Main plot - rainfall-runoff histograms ======================================
main_plot <- tidy_boxcox_streamflow |>
  ggplot(aes(x = rainfall, y = boxcox_streamflow, colour = control_or_change, fill = control_or_change)) +
  geom_point(shape = 21, alpha = 0.7) +
  geom_smooth(formula = y ~ x, method = lm, se = FALSE, linewidth = 0.25) +
  labs(
    x = "Total Annual Precipitation (mm)",
    y = "Total Annual Streamflow (Box-Cox Transformed)"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_manual(values = c("#f5a6a7", "#d4e5f2")) +
  theme_bw() +
  facet_wrap(
    ~parameter,
    nrow = 3, ncol = 2, axis.labels = "margins"
  ) +
  theme(
    legend.title = element_blank(),
    strip.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.08, 0.95),
    legend.background = element_rect(fill = NULL, colour = "black", linewidth = 0.2),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, linewidth = 0.5)))



abc_labels <- data.frame(
  label = paste0(letters[1:length(control_parameters)], ")"), 
  parameter = unique(tidy_boxcox_streamflow$parameter)) |>
  geom_text(
    mapping = aes(
      x = max(tidy_boxcox_streamflow$rainfall),
      y = max(tidy_boxcox_streamflow$boxcox_streamflow),
      label = label
    ),
    inherit.aes = FALSE,
    size = 6
  ) # , fontface = "bold"



rainfall_runoff_plot <- main_plot + abc_labels

# Save results -----------------------------------------------------------------
## Rainfall-runoff relationship ================================================
ggsave(paste0("./Graphs/rainfall_sensitivity_bc_rainfall_runoff_a4_page_", get_date(), ".pdf"),
       plot = rainfall_runoff_plot,
       device = "pdf",
       units = "mm",
       width = 210,
       height = 210
       ) 



# Comparing rainfall bc changes and streamflow bc changes ----------------------
rainfall_runoff_streamflow_changes <- read_csv("./Results/rainfall_runoff_streamflow_changes.csv", show_col_types = FALSE)

# I am only interested in comparing the mean - intercept and sd and sd
tidy_rainfall_runoff_rainfall_changes <- tidy_boxcox_streamflow |> 
  filter(parameter %in% c("mean", "standard_deviation")) |> 
  add_column(model = "Synthetic Streamflow Model", .before = 1) 


tidy_rainfall_runoff_streamflow_changes <- rainfall_runoff_streamflow_changes |>
  filter(parameter %in% c("intercept", "standard_deviation")) |> 
  add_column(model = "Stochastic Rainfall Model", .before = 1) |> 
  select(colnames(tidy_rainfall_runoff_rainfall_changes))


streamflow_comparison <- rbind(tidy_rainfall_runoff_streamflow_changes, tidy_rainfall_runoff_rainfall_changes) |> 
  mutate(
    change_type = if_else(parameter == "standard_deviation", "Change in Peaks and Troughs", "Change in Vertical Axis"),
    parameter = if_else((parameter == "standard_deviation") & (model == "Stochastic Rainfall Model"), "standard_deviation_r", parameter)
  )

# Plot -------------------------------------------------------------------------

y_label_locations <- streamflow_comparison |>
  summarise(
    max_streamflow = max(streamflow) * 0.95,
    .by = change_type
  ) |>
  pull(max_streamflow) |>
  rep(each = 2) |>
  rev()

streamflow_comparison_plot <- streamflow_comparison |>
  ggplot(aes(x = time, y = streamflow, colour = control_or_change)) +
  geom_line() +
  labs(
    x = "Year(s) of Data",
    y = "Total Annual Streamflow (mm)"
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  facet_grid(change_type ~ model, scales = "free_y") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 10)
  ) 




abc_labels_streamflow <- data.frame(
  label = paste0(letters[1:length(unique(streamflow_comparison$parameter))], ")"),
  change_type = c("Change in Peaks and Troughs", "Change in Peaks and Troughs", "Change in Vertical Axis", "Change in Vertical Axis"),
  model = c("Stochastic Rainfall Model", "Synthetic Streamflow Model", "Stochastic Rainfall Model", "Synthetic Streamflow Model"),
  x_location = min(streamflow_comparison$time),
  y_location = y_label_locations
) |>
  geom_text(
    mapping = aes(
      x = x_location,
      y = y_location,
      label = label
    ),
    inherit.aes = FALSE,
    size = 6
  ) 


final_streamflow_comparison_plot <- streamflow_comparison_plot + abc_labels_streamflow

ggsave(paste0("./Graphs/streamflow_comparison_", get_date(), ".pdf"),
       plot = final_streamflow_comparison_plot,
       device = "pdf",
       width = 210,
       height = 135,
       units = "mm"
       )








