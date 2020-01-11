#####################################################
# Implementation of Focused BH for outer nodes filter
#####################################################

# Run Focused BH
# P: vector of p-values
# G: graph object (see Test_FocusedBH)
# q: target FDR level
# gamma: numerator function for FDP-hat; 
#         equal to gamma(t) = m*t by default
FocusedBH = function(P, G, q, gamma = NULL){
  m = G$m
  if(is.null(gamma)){
    gamma = function(p)(m*p)
  }
  num_rejections = get_num_outer_nodes(P, G)
  P_aug = sort(c(0,P), decreasing = TRUE)
  FDP_hat = gamma(P_aug)/num_rejections
  FDP_hat[is.na(FDP_hat)] = 0
  FDP_hat[G$m+1] = 0
  idx = min(which(FDP_hat <= q))
  R = P <= P_aug[idx]
  return(R)
}

# Compute outer nodes for a set of nodes in a graph
# S: logical inclusion vector for set of nodes
# G: graph object
get_outer_nodes = function(S, G){
  m = G$m
  Pa = G$Pa
  marked = logical(m)
  not_outer = logical(m)
  state = c()
  state$marked = marked
  state$not_outer = not_outer
  traverse = function(node, outer, state){
    if(!outer){
      state$not_outer[node] = TRUE
    }
    state$marked[node] = TRUE
    for(parent in Pa[[node]]){
      if(!state$not_outer[parent]){
        state$not_outer[parent] = TRUE
      }
      if(!state$marked[parent]){
        state = traverse(parent, FALSE, state)
      }
    }
    return(state)
  }
  for(node in which(S)){
    state = traverse(node, TRUE, state)
  }
  return(S & !state$not_outer)
}

# Get number of outer nodes for each possible p-value threshold
# P: vector of p-values
# G: graph object
get_num_outer_nodes = function(P, G){
  m = G$m
  Pa = G$Pa
  C = G$C
  
  # initialize things
  num_rejections = integer(m+1)
  rejections = sum(sapply(C, length) == 0) 
  num_rejections[1] = rejections
  
  in_rejection_set = !logical(m)
  
  ord = order(P, decreasing = TRUE)
  
  tol = 0
  # iteration k means you remove kth largest p-value
  for(k in 1:m){
    P_to_remove = P[ord[k]]
    # keep removing nodes with equal p-values from rejection set
    # and rewiring the graph
    for(j in k:m){
      # find next candidate node to remove
      to_remove = ord[j]
      # stopping condition
      if(!in_rejection_set[to_remove] | P[to_remove] < P_to_remove - tol){
        break
      }
      # remove the node
      in_rejection_set[to_remove] = FALSE
      # check if we are removing an outer node
      if(length(C[[to_remove]]) == 0){
        rejections = rejections - 1
      }
      # get active parents and children of the node
      active_parents = Pa[[to_remove]][in_rejection_set[Pa[[to_remove]]]]
      active_children = C[[to_remove]][in_rejection_set[C[[to_remove]]]]
      # update lists of parents and children
      for(parent in active_parents){
        C[[parent]] = setdiff(C[[parent]], to_remove)
      }
      for(child in active_children){
        Pa[[child]] = setdiff(Pa[[child]], to_remove)
      }
      # rewire the graph
      for(parent in active_parents){
        for(child in active_children){
          Pa[[child]] = union(Pa[[child]], parent)
          C[[parent]] = union(C[[parent]], child)
        }
        # check if we are adding an outer node
        if(length(C[[parent]]) == 0){
          rejections = rejections + 1
        }
      }
    }
    num_rejections[k+1] = rejections
  }
  return(num_rejections)
}