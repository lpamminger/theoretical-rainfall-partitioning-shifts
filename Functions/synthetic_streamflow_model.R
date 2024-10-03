synthetic_streamflow_model <- function(control_parameters, control_rainfall, random) {
  
  force(control_parameters)
  force(control_rainfall)
  force(random)
  
  
  # set a constant seed if random is FALSE or generate random seed if true
  if (random == FALSE) {
    set.seed(8675309) 
  } else if (random == TRUE) {
    set.seed(NULL)
  }
  
  
  # get control_parameters
  control_intercept <- control_parameters[1]
  control_slope <- control_parameters[2]
  control_autocorr <- control_parameters[3]
  control_stand_dev <- control_parameters[4]
  control_skew <- control_parameters[5]
  
  # control decay value
  control_decay_value <- control_intercept + (control_slope * mean(control_rainfall))
  
  # set previous value
  control_boxcox_streamflow_previous <- control_decay_value
  
  control_boxcox_streamflow <- numeric(length = length(control_rainfall))
  
  # Do I care about the mean-boxcox streamflow? not including it may change some things
  
  # get standard deviation and skewness to the skewed normal parameters
  # The error does not change.
  control_skewed_normal_parameters <- to_skewed_normal_parameters(
                                        user_mean = 0,
                                        user_std = control_stand_dev,
                                        user_skewness = control_skew)
  # Run the model
  for (time_step in seq_along(control_rainfall)) {
    
    control_boxcox_streamflow[time_step] <- control_intercept + 
      (control_slope * control_rainfall[time_step]) + 
      (control_autocorr * (control_boxcox_streamflow_previous - control_decay_value)) +
      rsn(n = 1, xi = control_skewed_normal_parameters[1], omega = control_skewed_normal_parameters[2], alpha = control_skewed_normal_parameters[3])
    
    control_boxcox_streamflow_previous <- control_boxcox_streamflow[time_step]
    
  }
  
  
  
  function(change_parameters, change_rainfall) {
    
    #browser() # for debugging
    
    force(change_parameters)
    force(change_rainfall)
    
    
    # Reset the seed
    if (random == FALSE) {
      set.seed(8675309) 
    } else if (random == TRUE) {
      set.seed(NULL)
    }
    
    
    # get change_parameters
    change_intercept <- change_parameters[1]
    change_slope <- change_parameters[2]
    change_autocorr <- change_parameters[3]
    change_stand_dev <- change_parameters[4]
    change_skew <- change_parameters[5]
    
    # change decay value
    change_decay_value <- change_intercept + (change_slope * mean(change_rainfall))
    
    # set previous value
    change_boxcox_streamflow_previous <- change_decay_value
    
    change_boxcox_streamflow <- numeric(length = length(change_rainfall))
    
    # Do I care about the mean-boxcox streamflow? not including it may change some things
    
    # get standard deviation and skewness to the skewed normal parameters
    change_skewed_normal_parameters <- to_skewed_normal_parameters(
                                         user_mean = 0,
                                         user_std = change_stand_dev,
                                         user_skewness = change_skew
                                         )
    
    # Run the model
    for (time_step in seq_along(change_rainfall)) {
      
      change_boxcox_streamflow[time_step] <- change_intercept + 
        (change_slope * change_rainfall[time_step]) + 
        (change_autocorr * (change_boxcox_streamflow_previous - change_decay_value)) +
        rsn(n = 1, xi = change_skewed_normal_parameters[1], omega = change_skewed_normal_parameters[2], alpha = change_skewed_normal_parameters[3])
      
      change_boxcox_streamflow_previous <- change_boxcox_streamflow[time_step]
      
    }
    

    # I want the function to return the mean and actual boxcox streamflow (as a tibble) |control_mean|control_actual|change_mean|change_actual
    result <- cbind(control_rainfall, control_boxcox_streamflow, change_rainfall, change_boxcox_streamflow)
    return(result)
    
  }
}