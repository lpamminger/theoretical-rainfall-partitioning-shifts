# Streamflow sensitivity analysis

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")
par(mfrow = c(1,1))


# Import libraries--------------------------------------------------------------
pacman::p_load(sn, moments, tidyverse)
# purrr and ggplot are part of the tidyverse
# stochastic_rainfall_generator function requires moments package
# ggplot2 and ggpubr for plotting (ggpubr is an extension of ggplot)
# synthetic_streamflow_model function requires sn package for skewed normal distribution




# Import functions -------------------------------------------------------------
source("./Functions/box_cox_transforms.R")
source("./Functions/synthetic_streamflow_model.R") 
source("./Functions/adjusting_parameters.R")
source("./Functions/modified_stochastic_rainfall_generator.R")
source("./Functions/utility.R")




# Length of timeseries ---------------------------------------------------------
skip <- 2
user_time_series_length <- 100 # the length of rainfall determines the length of streamflow
time_series_length <- user_time_series_length + skip # we remove the first two values as the model spins up

# Generate rainfall using a annual stochastic rainfall model -------------------
generated_rainfall <- modified_stochastic_rainfall_generator(
                        parameter_vector = c(
                          "mean" = 1006, 
                          "sd" = 221, 
                          "auto" = 0.015, 
                          "skew" = 0.19), 
                        length_of_generated_rainfall = time_series_length,
                        random = FALSE)


# Generate box-cox streamflow using synthetic model ----------------------------

## Making up control parameters ================================================
control_parameters <- c(-4.1, 0.017, 0.16, 2, 0.014) # c(a0, a1, a2, a3, a4) c(3, 0.01, 0.12, 0.45, -0.9)
### Guidelines for making the control parameters:
### - typical values for a2 (autocorr.) is 0.05 and 0.1 for Australia
### - streamflow cannot be greater than precipiation
### - streamflow cannot be less than zero


## Make up change parameters ===================================================
multipliers_for_control_parameters <- c(1.5, 1.3, 3, 2, 65) 

change_parameters <- imap(.x = multipliers_for_control_parameters, 
                          .f = change_parameter_set_function, 
                          control_parameter_set = control_parameters)


## Combine control and change parameter sets in a single list ==================
parameter_list <- c(list(control_parameters), change_parameters)
names(parameter_list) <- c("control",
                           "intercept",
                           "slope",
                           "autocorrelation",
                           "standard_deviation",
                           "skewness")




## Apply model =================================================================

### Part 1 of function factory - store the control values
change_synthetic_streamflow_model <- synthetic_streamflow_model(
                                       control_parameters = parameter_list$control,
                                       control_rainfall = generated_rainfall,
                                       random = FALSE
                                       )

### Part 2 of function factory - add the change values
boxcox_streamflow_list <- map(.x = parameter_list[-1], # exclude control
                              .f = change_synthetic_streamflow_model, 
                              change_rainfall = generated_rainfall
                              )

# add another column - discard the first two columns as the assumed values have influence on slope and intercept
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
                         .before = 1)



control_minus_change_residual <- function(control_rainfall, control_boxcox_streamflow, change_boxcox_streamflow) {
  
  control_line_of_best_fit <- lm(control_boxcox_streamflow ~ control_rainfall)$fitted.values
  
  control_residual <- control_boxcox_streamflow - control_line_of_best_fit
  
  change_residual <- change_boxcox_streamflow - control_line_of_best_fit
  
  result <- list("control_residual" = control_residual,
                 "change_residual" = change_residual
  )
  
  return(result)
  
}

dot_to_line_residual <- function(rainfall, boxcox_streamflow) {
  line_of_best_fit <- lm(boxcox_streamflow ~ rainfall)$fitted.values
  residual <- boxcox_streamflow - line_of_best_fit
}


boxcox_streamflow <- boxcox_streamflow |> 
                       mutate(
                         control_residual = control_minus_change_residual(
                           control_rainfall = control_rainfall,
                           control_boxcox_streamflow = control_boxcox_streamflow,
                           change_boxcox_streamflow = change_boxcox_streamflow
                           )[[1]],
                         change_residual = control_minus_change_residual(
                           control_rainfall = control_rainfall,
                           control_boxcox_streamflow = control_boxcox_streamflow,
                           change_boxcox_streamflow = change_boxcox_streamflow
                         )[[2]],
                         control_dot_to_line = dot_to_line_residual(
                           rainfall = control_rainfall, 
                           boxcox_streamflow = control_boxcox_streamflow
                           ),
                         change_dot_to_line = dot_to_line_residual(
                           rainfall = change_rainfall, 
                           boxcox_streamflow = change_boxcox_streamflow
                           ),
                         .by = parameter
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


value_to_name <- c("rainfall", "boxcox_streamflow", "boxcox_streamflow_residual", "boxcox_dot_to_line_residual")

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
                              mean = mean(boxcox_dot_to_line_residual),
                              median = median(boxcox_dot_to_line_residual),
                              sd = sd(boxcox_dot_to_line_residual),
                              skew = skewness(boxcox_dot_to_line_residual),
                              auto = acf(boxcox_dot_to_line_residual, lag.max = 1, plot = FALSE)$acf[2],
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
    "intercept",
    "slope",
    "autocorrelation",
    "standard_deviation",
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
    legend.position.inside = c(0.75, 0.15),
    legend.background = element_rect(fill = NULL, colour = "black", linewidth = 0.2),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 14)
  ) + 
  guides(colour = guide_legend(override.aes = list(size = 5, linewidth = 0.5)))



## Inset plots - Residual histograms ===========================================
### annotate_custom does not allow multiple annotations for each panel

# calculate breaks for every histogram to make sure they all look the same
x <- tidy_boxcox_streamflow |>
  filter(parameter == "intercept") |>
  filter(control_or_change == "Control") |> 
  pull(boxcox_streamflow_residual) 

y <- tidy_boxcox_streamflow |>
  filter(parameter == "intercept") |>
  filter(control_or_change == "Change") |> 
  pull(boxcox_streamflow_residual) 

# bins should be offset by:
binwidth_offset <- abs(unique(y - x)[1] / 2)



# Add a different inset histogram for each panel
# Taken from https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
annotation_custom2 <- function(grob, data, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) {
  layer(
    data = data, stat = StatIdentity, position = PositionIdentity,
    geom = ggplot2:::GeomCustomAnn,
    inherit.aes = TRUE, params = list(
      grob = grob,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  )
}





histogram_generator <- function(parameter_specific_data) {
  parameter_specific_data |>
    ggplot(aes(x = boxcox_streamflow_residual, fill = control_or_change, colour = control_or_change)) +
    geom_histogram(position = "identity", binwidth = binwidth_offset, alpha = 0.4, show.legend = FALSE, linewidth = 0.1) +
    theme_bw() +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x = "Residual",
      y = "Frequencey"
    ) +
    theme(
      axis.title = element_text(size = 9),
      plot.background = element_blank(),
      axis.text = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
}





make_annotation_inset <- function(parameter, data) {
  
  # extract
  parameter_specific_data <- data |> filter(parameter == {{ parameter }})

  # make histogram and turn it into a grob
  grob_histogram <- histogram_generator(parameter_specific_data) |> 
    ggplotGrob()

  # convert grob into annotation
  annotation <- annotation_custom2(
    grob = grob_histogram, 
    data = parameter_specific_data, 
    xmin = 410,
    xmax = 1000,
    ymin = 19,
    ymax = 32
  )
  
  return(annotation)
}


# Inset histogram needs to be scaled correctly
# Binwidth should be constant. I will have to wrestle the bins. bins for intercept does not look correct




inset_histograms <- map(
  .x = unique(tidy_boxcox_streamflow$parameter),
  .f = make_annotation_inset,
  data = tidy_boxcox_streamflow
)

abc_labels <- data.frame(label = paste0(letters[1:length(control_parameters)],")"), parameter = unique(tidy_boxcox_streamflow$parameter)) |> 
  geom_text(
  mapping = aes(
    x = max(tidy_boxcox_streamflow$rainfall), 
    y = max(tidy_boxcox_streamflow$boxcox_streamflow), 
    label = label
    ), 
  inherit.aes = FALSE, 
  size = 6) #, fontface = "bold"


rainfall_runoff_plot <- main_plot + inset_histograms + abc_labels






# Lag-1 streamflow vs streamflow -----------------------------------------------
lag_1_streamflow_graph <- tidy_boxcox_streamflow |>
  filter(parameter == "autocorrelation") |>
  mutate(
    lag_boxcox_streamflow = lag(boxcox_streamflow),
    .by = control_or_change
  ) |>
  tail(-2) |> # remove first two values of tibble because of lag
  ggplot(aes(x = boxcox_streamflow, y = lag_boxcox_streamflow, colour = control_or_change, fill = control_or_change)) +
  geom_point(shape = 21, alpha = 0.7) +
  geom_smooth(formula = y ~ x, method = lm, se = FALSE, linewidth = 0.25) +
  labs(
    x = "Total Annual Streamflow (Box-Cox Transformed)",
    y = "Lag-1 Total Annual Streamflow (Box-Cox Transformed)"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_manual(values = c("#f5a6a7", "#d4e5f2")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.89),
    legend.background = element_rect(fill = NULL, colour = "black", linewidth = 0.2),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2, linewidth = 0.5)))







# Save results -----------------------------------------------------------------
write_csv(tidy_boxcox_streamflow, "./Results/rainfall_runoff_streamflow_changes.csv")

## Rainfall-runoff relationship ================================================
ggsave(paste0("./Graphs/streamflow_sensitivity_bc_rainfall_runoff_a4_page_", get_date(), ".pdf"),
       plot = rainfall_runoff_plot,
       device = cairo_pdf,
       units = "mm",
       width = 210,
       height = 297) 




## Lag-1-Streamflow ============================================================
ggsave(paste0("./Graphs/lag_1_streamflow_graph_", get_date(), ".pdf"),
  plot = lag_1_streamflow_graph,
  device = "pdf",
  width = 100,
  height = 100,
  units = "mm"
)



