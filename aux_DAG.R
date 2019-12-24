FocusedBH_outer_nodes = function(P, G, q, gamma = NULL){
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

FocusedBH = function(filter_name, P, G, q, gamma = NULL){
  cat(sprintf("Running Focused BH...\n"))
  if(filter_name == "outer_nodes"){
    return(FocusedBH_outer_nodes(P, G, q, gamma))
  }
  m = G$m
  if(is.null(gamma)){
    gamma = function(p)(m*p)
  }
  ord = order(P)
  fdr_thresh = max(P[p.adjust(P, "fdr") <= q])
  k_start = sum(P <= fdr_thresh)
  
  if(k_start == 0){
    return(logical(m))
  } else{
    for(k in k_start:1){
      cat(sprintf("Considering rejection set of size %d...\n", k))
      R = logical(m)
      R[ord[1:k]] = TRUE
      R_F = run_filter(P, R, G, filter_name)
      FDP_hat = gamma(P[ord[k]])/sum(R_F)
      if(FDP_hat <= q){
        return(R)
      }
    }
    return(logical(m))
  }
}

run_filter = function(P, R, G, filter_name){
  if(filter_name == "outer_nodes"){
    return(get_outer_nodes(R, G))
  }
  
  if(filter_name == "REVIGO"){
    R = R & P < 0.5
    goterms = tibble(names(which(R)), P[R]) %>% format_tsv(col_names = FALSE)
    command = sprintf("echo \"%s\" | java -jar REVIGO/RevigoStandalone.jar /dev/stdin --stdout --cutoff=0.5", goterms)
    output = system(command, intern = TRUE)
    matrix_output = do.call("rbind", lapply(output[3:length(output)], function(str)(unname(strsplit(str, split = "\t"))[[1]])))
    df_output = as.data.frame(matrix_output[2:nrow(matrix_output),,drop=F])
    names(df_output) = matrix_output[1,]
    R_F = R
    R_F[which(R_F)] = df_output$eliminated == 0
    
    return(R_F)
  }
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

get_ancestors = function(R, G){
  n = G$n
  Pa = G$Pa
  traverse = function(node, ancestor){
    if(!ancestor[node]){
      ancestor[node] = TRUE
      for(parent in Pa[[node]]){
        ancestor = traverse(parent, ancestor)
      }
    }
    return(ancestor)
  }
  ancestor = logical(n)
  for(node in which(R)){
    ancestor = traverse(node, ancestor)
  }
  return(ancestor)
}

get_descendants = function(R, G){
  n = G$n
  C = G$C
  traverse = function(node, descendant){
    if(!descendant[node]){
      descendant[node] = TRUE
      for(child in C[[node]]){
        descendant = traverse(child, descendant)
      }
    }
    return(descendant)
  }
  descendant = logical(n)
  for(node in which(R)){
    descendant = traverse(node, descendant)
  }
  return(descendant)
}

subset_graph = function(to_remove, G){
  C = G$C
  Pa = G$Pa
  for(i in to_remove){
    # print(i)
    for(j in C[[i]]){
      Pa[[j]] = union(setdiff(Pa[[j]], i), Pa[[i]]) # remove i, add i's parents
    }
    for(j in Pa[[i]]){
      C[[j]] = union(setdiff(C[[j]], i), C[[i]]) # remove i, add i's children
    }
    C[[i]] = NULL
    Pa[[i]] = NULL
  }
  G_new = c()
  G_new$C = C
  G_new$Pa = Pa
  return(G_new)
}

get_depths = function(Pa){
  m = length(Pa)
  depths = numeric(m)
  while(TRUE){
    done = TRUE
    for(j in 1:m){
      if(depths[j] == 0){
        done = FALSE
        if(length(Pa[[j]]) == 0){
          depths[j] = 1
        } 
        else{
          if(all(depths[Pa[[j]]] > 0)){
            depths[j] = 1 + max(depths[Pa[[j]]])
          }
        }
      }
    }
    if(done){
      break
    }
  }
  return(depths)  
}

get_heights = function(C){
  m = length(C)
  heights = numeric(m)
  while(TRUE){
    done = TRUE
    for(j in 1:m){
      if(heights[j] == 0){
        done = FALSE
        if(length(C[[j]]) == 0){
          heights[j] = 1
        } 
        else{
          if(all(heights[C[[j]]] > 0)){
            heights[j] = 1 + max(heights[C[[j]]])
          }
        }
      }
    }
    if(done){
      break
    }
  }
  return(heights)  
}
