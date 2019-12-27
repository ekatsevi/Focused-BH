tree_file = sprintf("%s/data/biobank/coding19.tsv", base_dir)
graph_file = sprintf("%s/data/processed/ICD_graph.Rda", base_dir)

# create graph structure
if(!file.exists(graph_file)){
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
  
  G = c()
  G$m = m
  G$C = C
  G$Pa = Pa
  G$depths = get_depths(Pa)
  
  # potentially filter G at this point
    
  leaves = sapply(G$C, length) == 0
  groups = vector("list", G$m)
  for(j in 1:G$m){
    R_delta = logical(G$m)
    R_delta[j] = TRUE
    groups[[j]] = which(get_descendants(R_delta, G) & leaves)
  }
  G$groups = groups
  
  save(G, file = graph_file)
}