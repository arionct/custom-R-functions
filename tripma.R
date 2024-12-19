tripma = function(Y, w, type='uniform', a=0.5) {
  mov = c()
  
  # no window for first w points
  for(k in 1:w) {
    mov = c(mov, NA)
  }
  
  if(type == 'uniform') {
    # uniform weights
    weights = rep(1/(w+1), w+1)
    
    for(k in (w+1):length(Y)) {
      mov = c(mov, mean(Y[(k-w):k]))
    }
  }
  
  else if(type == 'exponential') {
    # exponential weights
    weights = a^(w:0)
    weights = weights / sum(weights)
    
    for(k in (w+1):length(Y)) {
      window = Y[(k-w):k]
      mov = c(mov, sum(weights * window))
    }
  }
  
  else {
    stop("Invalid 'type' argument. Choose 'uniform' or 'exponential'.")
  }
  
  # objects to return
  results = list('approx'  = mov,
                 'weights' = weights)
  
  return(results)
}