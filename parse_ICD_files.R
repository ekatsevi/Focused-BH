tree_file = sprintf("%s/data/biobank/coding19.tsv", base_dir)
graph_file = sprintf("%s/data/processed/ICD_graph.Rda", base_dir)

# create graph structure
if(!file.exists(graph_file)){
  cat(sprintf("Parsing ICD tree...\n"))
  tree = read_tsv(tree_file, col_types = c("cciic"))
  m = nrow(tree)
  id2idx = list()
  for(j in 1:m){
    id2idx[as.character(tree$node_id[j])] = j
  }
  C = vector("list", m)
  Pa = vector("list", m)
  for(j in 1:m){
    node_id = tree$node_id[j]
    parent_id = tree$parent_id[j]
    if(parent_id != 0){
      parent_idx = id2idx[[as.character(parent_id)]]
      Pa[[j]] = c(Pa[[j]], parent_idx)
      C[[parent_idx]] = c(C[[parent_idx]], j)
    }
  }
  
  node_names = tree$coding

  leaves = sapply(C, length) == 0
  item_names = node_names[leaves]
  num_items = sum(leaves)
  
  item_index = integer(G$num_items)
  item_index[leaves] = 1:num_items
  sets = vector("list", m)
  for(node in 1:m){
    R_delta = logical(m)
    R_delta[node] = TRUE
    sets[[node]] = item_index[get_descendants(R_delta, C) & leaves]
  }
  
  G = c()
  G$m = m
  G$num_items = num_items
  G$C = C
  G$Pa = Pa
  G$sets = sets
  G$node_names = node_names
  G$item_names = item_names
  
  # convert graph to Structured Holm format
  G_SH = new("DAGstructure", parents = G$Pa, children = G$C, sets = G$sets, twoway = FALSE)
  
  # convert graph structure to Yekutieli format
  roots = which(sapply(G$Pa, length) == 0)
  # add an imaginary root node, whose roots are current root nodes
  G_augmented = G
  G_augmented$m = G_augmented$m + 1
  G_augmented$C[[m+1]] = roots
  G_augmented$Pa[[m+1]] = integer(0)
  for(root in roots){
    G_augmented$Pa[[root]] = m + 1
  }
  
  G_Yekutieli = t(do.call(cbind, 
                          lapply(1:G_augmented$m, 
                                 function(i)(if(length(G_augmented$C[[i]]) > 0) sapply(G_augmented$C[[i]], 
                                                                                     function(j)(c(as.character(i), as.character(j)))) else matrix(0,2,0)))))
  save(G, G_SH, G_Yekutieli, file = graph_file)
}