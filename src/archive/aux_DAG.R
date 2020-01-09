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

subset_graph = function(subset, G){
  n_new = sum(subset)
  children_new = c()
  groups_new = c()
  node_indices = cumsum(subset)
  i_new = 1
  for(i in 1:G$n){
    if(subset[i]){
      children = G$C[[i]]
      group = G$groups[[i]]
      children = node_indices[children[subset[children]]]
      children_new[[i_new]] = children
      groups_new[[i_new]] = group
      i_new = i_new + 1
    }
  }
  # get parents
  parents_new = vector("list", n_new)
  for(i_new in 1:n_new){
    for(j_new in children_new[[i_new]]){
      parents_new[[j_new]] = c(parents_new[[j_new]], i_new)
    }
  }
  G_new = c()
  G_new$n = n_new
  G_new$groups = groups_new
  G_new$C = children_new
  G_new$Pa = parents_new
  return(G_new)
}

subset_graph_2 = function(to_remove, G){
  C = G$C
  Pa = G$Pa
  for(i in to_remove){
    print(i)
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

# # temporary
# REVIGO = function(P, R){
#   if(filter_name == "REVIGO"){
#     R = R & P < 0.5
#     if(sum(R) == 0){
#       return(NULL)
#     } else{
#       # goterms = tibble(names(which(R)), P[R]) %>% format_tsv(col_names = FALSE)
#       goterms = tibble(names(which(R))) %>% format_tsv(col_names = FALSE)
#       command = sprintf("echo \"%s\" | java -jar REVIGO/RevigoStandalone.jar /dev/stdin --stdout --cutoff=0.5", goterms)
#       output = system(command, intern = TRUE)
#       matrix_output = do.call("rbind", 
#                               lapply(output[3:length(output)], 
#                                      function(str)(unname(strsplit(str, split = "\t"))[[1]])))
#       df_output = as.data.frame(matrix_output[2:nrow(matrix_output),,drop=F])
#       names(df_output) = matrix_output[1,]
#       df_output$frequency = sapply(as.character(df_output$frequency), function(freq)(strsplit(freq, split = "%")[[1]]))
#       return(df_output %>% as_tibble() %>%
#                select("term_ID", "frequency", "log10_p-value", 
#                       "uniqueness", "dispensability", "representative", "eliminated") %>%
#                mutate_all(as.character) %>% 
#                mutate_at(c("frequency", "log10_p-value", "uniqueness", 
#                            "dispensability", "eliminated"),
#                          as.numeric))
#       df_output = df_output %>% as_tibble() %>%
#         dplyr::select("term_ID", "frequency", 
#                       "uniqueness", "dispensability", "representative", "eliminated") %>%
#         mutate_all(as.character) %>% 
#         mutate_at(c("frequency", "uniqueness", 
#                     "dispensability", "eliminated"),
#                   as.numeric)
#     }
#   }
# }
