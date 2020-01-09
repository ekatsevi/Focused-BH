# Focused BH 
# filter_name: name of the filter to be applied
# P:           vector of m p-values
# G:           graph structure
# q:           target FDR level
# V_hats:      list of functions to estimate number of false discoveries
FocusedBH = function(filter_name, P, G, q, V_hats = NULL){
  # apply accelerated function for special case of outer nodes filter
  if(filter_name == "outer_nodes"){
    return(FocusedBH_outer_nodes(P, G, q, V_hats))
  }
  cat(sprintf("Running Focused BH...\n"))
  
  # parse V_hats based on how many were given (0, 1, or more than 1)
  m = G$m
  V_hats_input = V_hats
  if(is.null(V_hats_input)){
    V_hats = list(function(t)(m*t))
    names(V_hats) = "Focused BH original"
  } else if(!is.list(V_hats_input)){
    V_hats = list()
    V_hats[[1]] = V_hats_input
    names(V_hats) = "Focused BH"
  } else{
    V_hats = V_hats_input
  }
  
  # set up boolean matrix for rejections for each V_hat
  num_V_hats = length(V_hats)
  rejections = matrix(FALSE, num_V_hats, m)
  
  # find BH threshold to start search for Focused BH
  BH_rejections = p.adjust(P, "fdr") <= q
  if(sum(BH_rejections) == 0){
    return(rejections)
  }
  BH_thresh = max(P[BH_rejections])
  
  # find number of rejections to start at
  k_start = sum(P <= BH_thresh)
  if(filter_name == "REVIGO"){
    k_start = min(k_start, 350)
  }

  # iteratively remove least significant hypothesis
  # and check if Focused BH estimate falls below q
  ord = order(P)
  complete = logical(num_V_hats) # see if all thresholds have been found
  for(k in k_start:1){
    cat(sprintf("Considering rejection set of size %d...\n", k))
    # define rejection set comprised of k most significant hypotheses
    R = logical(m)
    R[ord[1:k]] = TRUE
    # apply filter to this rejection set
    R_F = run_filter(P, R, G, filter_name)
    # iterate over definitions of V_hat
    for(index in which(!complete)){
      # check if FDP_hat is below the threshold
      FDP_hat = V_hats[[index]](P[ord[k]])/sum(R_F)
      if(FDP_hat <= q){
        rejections[index,] = R
        cat(sprintf("Found threshold for %s!\n", names(V_hats)[index]))
        complete[index] = TRUE
        if(all(complete)){
          if(num_V_hats == 1){
            rejections = as.logical(rejections)
          }
          return(rejections)
        }
      }
    }
  }
  if(num_V_hats == 1){
    rejections = as.logical(rejections)
  }
  return(rejections)
}

# Focused BH, accelerated for outer nodes filter
# P:           vector of m p-values
# G:           graph structure
# q:           target FDR level
# V_hats:      list of functions to estimate number of false discoveries
FocusedBH_outer_nodes = function(P, G, q, V_hats = NULL){
  # parse V_hats based on how many were given (0, 1, or more than 1)
  m = G$m
  V_hats_input = V_hats
  if(is.null(V_hats_input)){
    V_hats = list(function(t)(m*t))
    names(V_hats) = "Focused BH original"
  } else if(!is.list(V_hats_input)){
    V_hats = list()
    V_hats[[1]] = V_hats_input
    names(V_hats) = "Focused BH"
  } else{
    V_hats = V_hats_input
  }
  num_V_hats = length(V_hats)
  
  # find BH threshold to start search for Focused BH
  rejections = matrix(FALSE, num_V_hats, m)
  BH_rejections = p.adjust(P, "fdr") <= q
  if(sum(BH_rejections) == 0){
    return(rejections)
  }
  BH_thresh = max(P[BH_rejections])
  
  # get number of outer nodes for every possible threshold
  num_rejections = get_num_outer_nodes(P, G)
  
  # find threshold for each definition of V_hat
  P_aug = sort(c(0, P[P <= BH_thresh]), decreasing = TRUE)
  for(index in 1:num_V_hats){
    V_hat = sapply(P_aug, V_hats[[index]])    
    FDP_hat = V_hat/num_rejections[(1 + (m+1) - length(P_aug)):(m+1)]
    FDP_hat[is.na(FDP_hat)] = 0
    idx = min(which(FDP_hat <= q))
    rejections[index,] = P <= P_aug[idx]
  }
  if(num_V_hats == 1){
    rejections = as.logical(rejections)
  }
  return(rejections)
}