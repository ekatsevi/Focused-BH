# make a line
make_DAG_line = function(D, root_node_size, leaf_node_size){
  node_sizes = round(seq(leaf_node_size, root_node_size, length.out = D))
  groups = sapply(1:D, function(j)(1:node_sizes[j]))
  DAG_info = c()
  DAG_info$n = D
  DAG_info$m = root_node_size
  DAG_info$groups = groups
  return(DAG_info)
}

# make a tree
make_DAG_tree = function(D, num_roots, num_children, leaf_node_size){
  if(length(num_children) == 1){
    num_children = rep(num_children, D-1)
  }
  layer_sizes = num_roots*cumprod(c(1,num_children))
  group_sizes = rev(leaf_node_size*cumprod(rev(c(num_children,1))))
  n = sum(layer_sizes)  
  m = num_roots*group_sizes[1]
  node_idx = 1
  groups = list()
  for(depth in 1:D){
    for(i in 1:layer_sizes[depth]){
      groups[[node_idx]] = ((i-1)/layer_sizes[depth]*m+1):(i/layer_sizes[depth]*m)
      node_idx = node_idx + 1
    }
  }
  DAG_info = c()
  DAG_info$n = n
  DAG_info$m = m
  DAG_info$groups = groups
  return(DAG_info)
}

make_DAG_forest_of_lines = function(D, W, root_node_size, leaf_node_size){
  node_sizes = round(seq(leaf_node_size, root_node_size, length.out = D))
  groups = list()
  for(w in 1:W){
    groups = c(groups, lapply(1:D, function(j)((w-1)*root_node_size + (1:node_sizes[j]))))
  }
  DAG_info = c()
  DAG_info$n = D*W
  DAG_info$m = root_node_size*W
  DAG_info$groups = groups
  return(DAG_info)
}

# change DAG from Jason's format to mine
make_DAG_GO = function(groups, children){
  genes = sort(unique(unlist(groups)))
  m = length(genes)
  gene_to_index = numeric(max(genes))
  for(i in 1:m){
    gene_to_index[genes[i]] = i
  }
  # note: removing duplicate groups here
  groups_new = unique(sapply(groups, function(group)(gene_to_index[group])))
  n = length(groups_new)
  # Adj = matrix(FALSE, n,n)
  # for(i in 1:n){
  #   Adj[i,setdiff(children[[i]], i)] = TRUE
  # }
  DAG_info = c()
  DAG_info$n = n
  DAG_info$m = m
  DAG_info$groups = groups_new
#  DAG_info$Adj = Adj
  return(DAG_info)
}

