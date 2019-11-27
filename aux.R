########################################
# Auxiliary functions for Test_FocusedBH
########################################

# Get the set of ancestors of a set of nodes 
# R: logical node inclusion vector
# G: graph object
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

# Get the set of descendants of a set of nodes 
# R: logical node inclusion vector
# G: graph object
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

# Construct a graph object for a subset of nodes
# subset: subset of nodes to restrict to
# G: graph object
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

# Obtain the depths of each node 
# Pa: list of parents
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

# Wrapper of igraph package to plot trees
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

# Plot rejections made by Focused BH in Test_FocusedBH
plot_rejections = function(R, G){
  ancestors = get_ancestors(R, G)
  nodes_to_plot = get_descendants(ancestors, G)
  G_to_plot = subset_graph(nodes_to_plot,G)
  vertex_shapes = "circle"
  frame_colors = NULL
  fill_colors = rep("white", G_to_plot$m)
  fill_colors[R[nodes_to_plot]] = "gray50"
  fill_colors[outer_nodes[nodes_to_plot]] = "red"
  vertex_labels = which(nodes_to_plot)
  vertex_size = 5
  plot_tree(G_to_plot, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)
}
