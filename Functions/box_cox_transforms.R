BCTransform <- function(y, lambda = 0, lambda_2 = 0) {
  if (lambda == 0L) { 
    log(y + lambda_2) 
  }
  else { 
    ((y + lambda_2)^lambda - 1) / lambda 
  }
}

BCTransformInverse <- function(yt, lambda = 0, lambda_2 = 0) {
  if (lambda == 0L) { 
    exp(yt) - lambda_2
  }
  
  else { 
    ((yt * lambda) + 1)^(1 / lambda) - lambda_2
  }
}