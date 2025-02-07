boxcox_transform <- function(y, lambda = 0, lambda_2 = 0) {
  # requires tidyverse

  # if(any(y + lambda_2 <= 0)){stop("Cannot transform value less than zero")}

  if (dplyr::near(lambda, rep(0, times = length(lambda)))) {
    log(y + lambda_2)
  } else {
    ((y + lambda_2)^lambda - 1) / lambda
  }
}




boxcox_inverse_transform <- function(yt, lambda = 0, lambda_2 = 0) {
  if (all(lambda == 0)) {
    exp(yt) - lambda_2
  } else {
    ((yt * lambda) + 1)^(1 / lambda) - lambda_2
  }
}



boxcox_lambda_generator <- function(precipitation, streamflow, lambda_2) {
  # , confidence_interval_percentage = 0.95
  ## Response variable q must be positive. Zero is not positive. Add a really small number + 1e-7.
  # Use MASS boxcox function to find best lambda
  boxcox_result_object <- MASS::boxcox(lm((streamflow + lambda_2) ~ precipitation, y = TRUE), # , y = TRUE IS NEEDED FOR WRAPPER TO WORK https://community.rstudio.com/t/error-in-eval-predvars-data-env-object-x-not-found-when-creating-a-function/129475/2
    lambda = seq(0, 5, 1 / 1000),
    plotit = FALSE
  )


  best_lambda <- boxcox_result_object$x[which.max(boxcox_result_object$y)]

  # Calculate confidence intervals
  ## What this code is doing
  ## Find the max likelihood. Then find the two 0.95 values. Return the lower and larger lambda values
  # confidence_interval <- range(boxcox_result_object$x[boxcox_result_object$y > max(boxcox_result_object$y) - qchisq(confidence_interval_percentage, df = 1) / 2])

  # Return c(best_lambda, lower_confidence_interval and upper_confidence_interval
  return(best_lambda)
}