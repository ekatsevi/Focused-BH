if(DAG_version == 1){
    children_filename = sprintf("%s/meta_gochild_dag.json", data_dir)
    groups_filename = sprintf("%s/meta_gogene_map.json", data_dir)
    children = read_json(children_filename, simplifyVector = TRUE)
    groups = read_json(groups_filename, simplifyVector = TRUE)
    genes = sort(unique(unlist(groups)))
    m = length(genes)
    gene_to_index = numeric(max(genes))
    for(i in 1:m){
      gene_to_index[genes[i]] = i
    }
    groups = sapply(groups, function(group)(gene_to_index[group]))  
    children = sapply(children, function(v)(v+1))
}

if(DAG_version == 2){
  D = 5
  num_roots = 1
  num_children = c(3,2,2,2)
  leaf_node_size = 1
  DAG_info = make_DAG_tree(D, num_roots, num_children, leaf_node_size)
  n = DAG_info$n
  m = DAG_info$m
  groups = DAG_info$groups
  children = NULL
}

if(DAG_version == 3){
  D = 4
  num_roots = 1
  num_children = 2
  leaf_node_size = 1
  DAG_info = make_DAG_tree(D, num_roots, num_children, leaf_node_size)
  n = DAG_info$n
  m = DAG_info$m
  groups = DAG_info$groups
  children = NULL
}

if(DAG_version == 4){
  D = 5
  num_roots = 20
  num_children = c(2,2,2,2)
  leaf_node_size = 1
  DAG_info = make_DAG_tree(D, num_roots, num_children, leaf_node_size)
  n = DAG_info$n
  m = DAG_info$m
  groups = DAG_info$groups
  children = NULL
}

if(DAG_version == 5){
  D = 6
  num_roots = 20
  num_children = c(2,2,2,2,2)
  leaf_node_size = 1
  DAG_info = make_DAG_tree(D, num_roots, num_children, leaf_node_size)
  n = DAG_info$n
  m = DAG_info$m
  groups = DAG_info$groups
  children = NULL
}


##################################
# G = make_DAG(groups, children)
# fill_colors = rep("white", G$n)
# fill_colors[!nulls_node] = "red"
# plot_tree_2(G, NULL, NULL, fill_colors, "", 5)

