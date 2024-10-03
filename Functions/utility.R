# Copying histogram bin function factory from advanced R 2e (put in separate folder)
binwidth_bins <- function(n){
  force(n)
  
  function(x){
    (max(x) - min(x)) / n
  }
}


# get date
get_date <- function() {
  str_remove_all(Sys.Date(), "-")
}


# Tom's start end splice
toms_start_end_splice <- function(input_vector) {
  
  NA_is_true <- !is.na(input_vector) # logicals of NA
  
  # Coherses logical to numeric. Makes start of NA -1 and end 1
  start_end_marking <- NA_is_true[2:NROW(NA_is_true)] - NA_is_true[1:(NROW(NA_is_true) - 1)] 
  
  if(NA_is_true[1] == TRUE) {
    index_start <- c(wf(NA_is_true, TRUE), (which(start_end_marking == 1) + 1)) # wf is which_first
  } else {
    index_start <- which(start_end_marking == 1) + 1
  }
  
  if(NA_is_true[NROW(NA_is_true)] == TRUE) {
    index_end <- c(which(start_end_marking == -1), wl(NA_is_true, TRUE)) # wl is which_last
    
  } else {
    index_end <- which(start_end_marking == -1)
  }
  
  start_end_table <- cbind(index_start,index_end)
  
} 

# Round any
round_any = function(x, accuracy, f = round){
  f(x / accuracy) * accuracy
}