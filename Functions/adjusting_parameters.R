# FUNCTION: to_skewed_normal_parameters
# 
# 
# Turning user defined skewness, variance and mean
# into parameters for skewed normal distribution. Equations found on wiki page: 
# https://en.wikipedia.org/wiki/Skew_normal_distribution

# INPUTS
# - user_skewness: a skewness between between -0.995 and 0.995. Limitation of equation.
# - user_std: user standard deviation
# - user_mean: user mean
#
#
# OUTPUTS
# - list of parameters for a skewed normal distribution  c(xi, omega, alpha)


to_skewed_normal_parameters <- function(user_mean, user_std, user_skewness){
  
  ## Step 1: Turn user defined skewness into alpha (shape parameter) and delta parameter
  
  # Define user skewness
  # Limitation: skewness must be between -0.995 and 0.995
  
  # Finding delta using methods of moments
  ## See estimation in wikipedia: https://en.wikipedia.org/wiki/Skew_normal_distribution
  abs_delta <- sqrt(((pi / 2) * abs(user_skewness)^(2 / 3)) / ((abs(user_skewness)^(2 / 3)) + (((4 - pi)/2)^(2 / 3))))
  delta <- sign(user_skewness) * abs_delta
  
  # Solving for alpha
  alpha <- delta / sqrt(1 - delta^2)
  
  ## Step 2: Turn user defined variance into omega (scale parameter)
  omega <- user_std / sqrt(1 - ((2 * delta^2) / pi))
  
  ## Step 3: Turn user defined mean into xi (location parameter)
  xi <- user_mean - (omega * delta * sqrt(2 / pi))
  
  return(c(xi, omega, alpha))
}








# DATE: 2023-05-10
# AUTHOR: Lucas Pamminger, Tim Peterson and Murray Peel
# VERSION: 1
#
# FUNCTION: synthetic_streamflow_model
# 
# Change a single value in the control parameters

# INPUTS
# - control_parameter_set: vector of control parameters
# - parameter_to_change: the position (index) of parameter in the control_parameter_set we want to change
# - multiplier: the multiplier applied to the specified parameter
# - 
#
# OUTPUTS
# - adjusted_parameter_set: a vector with adjusted parameters

change_parameter_set_function <- function(multiplier, parameter_to_change, control_parameter_set){
  
  adjusted_parameter_set <- control_parameter_set
  adjusted_parameter_set[parameter_to_change] <- adjusted_parameter_set[parameter_to_change] * multiplier
  
  return(adjusted_parameter_set)  
}