FocusedBH_outer_nodes = function(P, G, q, V_hats = NULL){
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

  rejections = matrix(FALSE, num_V_hats, m)
  
  BH_rejections = p.adjust(P, "fdr") <= q
  if(sum(BH_rejections) == 0){
    return(rejections)
  }
  
  BH_thresh = max(P[BH_rejections])
  
  num_rejections = get_num_outer_nodes(P, G)
  P_aug = sort(c(0, P[P <= BH_thresh]), decreasing = TRUE)
  
  for(index in 1:num_V_hats){
    V_hat = sapply(P_aug, V_hats[[index]])    
    FDP_hat = V_hat/num_rejections[(1 + (m+1) - length(P_aug)):(m+1)]
    FDP_hat[is.na(FDP_hat)] = 0
    idx = min(which(FDP_hat <= q))
    rejections[index,] = P <= P_aug[idx]
  }
  if(num_V_hats == 1){
    rejections = as.numeric(rejections)
  }
  return(rejections)
}

FocusedBH = function(filter_name, P, G, q, V_hats = NULL){
  if(filter_name == "outer_nodes"){
    return(FocusedBH_outer_nodes(P, G, q, V_hats))
  }
  cat(sprintf("Running Focused BH...\n"))
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
  rejections = matrix(FALSE, num_V_hats, m)
  
  BH_rejections = p.adjust(P, "fdr") <= q
  if(sum(BH_rejections) == 0){
    return(rejections)
  }
  
  BH_thresh = max(P[BH_rejections])
  
  ord = order(P)
  k_start = sum(P <= BH_thresh)
  
  if(filter_name == "REVIGO"){
    k_start = min(k_start, 350)
  }
  
  complete = logical(num_V_hats)
  for(k in k_start:1){
    cat(sprintf("Considering rejection set of size %d...\n", k))
    R = logical(m)
    R[ord[1:k]] = TRUE
    R_F = run_filter(P, R, G, filter_name)
    for(index in which(!complete)){
      FDP_hat = V_hats[[index]](P[ord[k]])/sum(R_F)
      if(FDP_hat <= q){
        rejections[index,] = R
        cat(sprintf("Found threshold for %s!\n", names(V_hats)[index]))
        complete[index] = TRUE
        if(all(complete)){
          if(num_V_hats == 1){
            rejections = as.numeric(rejections)
          }
          return(rejections)
        }
      }
    }
  }
  if(num_V_hats == 1){
    rejections = as.numeric(rejections)
  }
  return(rejections)
}

run_filter = function(P, R, G, filter_name){
  if(filter_name == "outer_nodes"){
    return(get_outer_nodes(R, G))
  }
  
  if(filter_name %in% c("REVIGO", "REVIGO2")){
    if(sum(R) == 0){
      R_F = R
    } else{
      if(filter_name == "REVIGO"){
        R = R & P < 0.5
        goterms = tibble(G$node_names[R], P[R]) %>% format_tsv(col_names = FALSE)        
      }
      if(filter_name == "REVIGO2"){
        goterms = tibble(G$node_names[R]) %>% format_tsv(col_names = FALSE)        
      }
      command = sprintf("echo \"%s\" | java -jar REVIGO/RevigoStandalone.jar /dev/stdin --stdout --cutoff=0.5", goterms)
      output = system(command, intern = TRUE)
      matrix_output = do.call("rbind", 
                              lapply(output[3:length(output)], 
                                     function(str)(unname(strsplit(str, split = "\t"))[[1]])))
      df_output = as.data.frame(matrix_output[2:nrow(matrix_output),,drop=F])
      names(df_output) = matrix_output[1,]
      df_output = df_output %>% 
        as_tibble() %>% 
        mutate(term_ID = as.character(term_ID), 
               eliminated = as.numeric(as.character(eliminated)))
      filtered_rejections = df_output %>% filter(eliminated == 0) %>% pull(term_ID)
      R_F = logical(length(R))
      R_F[G$node_names %in% filtered_rejections] = TRUE
    }
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

get_descendants = function(R, C){
  m = length(C)
  traverse = function(node, descendant){
    if(!descendant[node]){
      descendant[node] = TRUE
      for(child in C[[node]]){
        descendant = traverse(child, descendant)
      }
    }
    return(descendant)
  }
  descendant = logical(m)
  for(node in which(R)){
    descendant = traverse(node, descendant)
  }
  return(descendant)
}

subset_graph = function(to_remove, G){
  C = G$C
  Pa = G$Pa
  
  # remove nodes one at a time, rewiring the graph 
  # to preserve DAG relationships 
  for(i in which(to_remove)){
    for(j in C[[i]]){
      Pa[[j]] = union(setdiff(Pa[[j]], i), Pa[[i]]) # remove i, add i's parents
    }
    for(j in Pa[[i]]){
      C[[j]] = union(setdiff(C[[j]], i), C[[i]]) # remove i, add i's children
    }
  }
  
  # remove nodes and re-index the graph
  m_new = sum(!to_remove) # number of remaining nodes
  index = integer(G$m)
  index[which(!to_remove)] = 1:m_new
  Pa_new = sapply(Pa[!to_remove], function(parents_list)(index[parents_list]))
  C_new = sapply(C[!to_remove], function(children_list)(index[children_list]))

  # extract subsets of node names and sets
  node_names_new = G$node_names[!to_remove]
  sets_new = G$sets[!to_remove]
  item_names_new = G$item_names
  
  # if necessary, re-index the item names
  num_items_new = length(unique(unlist(sets_new)))
  if(num_items_new < G$num_items){
    remaining_items = unique(unlist(sets_new))
    item_names_new = G$item_names[remaining_items]
    item_index = integer(G$num_items)
    item_index[remaining_items] = 1:num_items_new
    sets_new = sapply(1:m_new, function(node)(item_index[sets_new[[node]]]))
  }
    
  # package into new graph object and return
  G_new = c()
  G_new$m = m_new
  G_new$num_items = num_items_new
  G_new$C = C_new
  G_new$Pa = Pa_new
  G_new$sets = sets_new
  G_new$node_names = node_names_new
  G_new$item_names = item_names_new
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

# get the number of descendants of each node that are
# among a list of marked nodes (only works for trees)
get_num_marked_descendants = function(C, marked){
  m = length(C)
  num_marked_descendants = rep(-1,m)
  while(TRUE){
    done = TRUE
    for(j in 1:m){
      if(num_marked_descendants[j] == -1){
        done = FALSE
        if(length(C[[j]]) == 0){
          num_marked_descendants[j] = as.numeric(marked[j])
        } 
        else{
          if(all(num_marked_descendants[C[[j]]] > -1)){
            num_marked_descendants[j] = as.numeric(marked[j]) + sum(num_marked_descendants[C[[j]]])
          }
        }
      }
    }
    if(done){
      break
    }
  }
  return(num_marked_descendants)  
}
