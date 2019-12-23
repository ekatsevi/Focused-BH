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
  m = G$m
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
  ancestor = logical(m)
  for(node in which(R)){
    ancestor = traverse(node, ancestor)
  }
  return(ancestor)
}

get_descendants = function(R, G){
  m = G$m
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
  descendant = logical(m)
  for(node in which(R)){
    descendant = traverse(node, descendant)
  }
  return(descendant)
}

subset_graph = function(subset, G){
  m_new = sum(subset)
  children_new = vector("list", m_new)
  node_indices = cumsum(subset)
  i_new = 1
  for(i in 1:G$m){
    if(subset[i]){
      children = G$C[[i]]
      children = node_indices[children[subset[children]]]
      children_new[[i_new]] = children
      i_new = i_new + 1
    }
  }
  # get parents
  parents_new = vector("list", m_new)
  for(i_new in 1:m_new){
    for(j_new in children_new[[i_new]]){
      parents_new[[j_new]] = c(parents_new[[j_new]], i_new)
    }
  }
  G_new = c()
  G_new$m = m_new
  G_new$C = children_new
  G_new$Pa = parents_new
  G_new$depths = get_depths(parents_new)
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
          if(depths[Pa[[j]]] > 0){
            depths[j] = 1 + depths[Pa[[j]]]
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

plot_tree = function(G, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size){
  library(igraph)
  m = G$m
  children = G$C
  parents = G$Pa
  edges = c()
  for(j in 1:m){
    for(child in children[[j]]){
      edges = c(edges, j, child)
    }
  }
  g<-graph(edges, n=m, directed=TRUE)
  roots = which(sapply(G$Pa, length) == 0)
  
  par(mar=c(0,0,1,0)) # make smaller margins
  plot(g, size = vertex_size, vertex.label = vertex_labels, vertex.size = vertex_size, vertex.color = fill_colors, vertex.shape = vertex_shapes,
       vertex.frame.color = frame_colors, layout = layout_as_tree(g, root = roots), asp = 0)
  par(mar=c(5.1, 4.1, 4.1, 2.1)) # reset margins
}