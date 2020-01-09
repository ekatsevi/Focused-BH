make_beta_random_locations = function(m, k){
  library(R.utils)
  beta = matrix(0, m, 1)
  nonzeros_tmp = c()
  if(k > 0){
    idx = withSeed(sample(1:m, k), seed = 572)
    beta[idx] = 1
  }
  return(beta)
}

make_beta_line = function(D, root_node_size, leaf_node_size, non_nulls_per_depth){
  m = root_node_size
  beta = matrix(0, m, 1)
  
  node_sizes = round(seq(leaf_node_size, root_node_size, length.out = D))
  groups = sapply(1:D, function(j)(1:node_sizes[j]))
  
  for(depth in 1:D){
    if(non_nulls_per_depth[depth] > 0){
      if(depth == D){
        new_nodes_at_depth = groups[[D-depth+1]]
      } else{
        new_nodes_at_depth = setdiff(groups[[D-depth+1]], groups[[D-depth]])
      }
      stopifnot(length(new_nodes_at_depth) >= non_nulls_per_depth[depth])
      beta[new_nodes_at_depth[1:non_nulls_per_depth[depth]]] = 1
    }
  }
  
  return(beta)
}

make_beta_forest_of_lines = function(D, W, root_node_size, leaf_node_size, first_non_null, num_non_nulls){
  m = W*root_node_size
  beta = matrix(0, m, 1)
  indices = 1:m
  beta[(indices-1) %% root_node_size + 1 >= first_non_null & (indices-1)/root_node_size < num_non_nulls] = 1
  return(beta)
}