suppressPackageStartupMessages(library(tidyverse))

### set top-level parameters
reps_per_experiment = 100 
reps = 100                         # Number of outer-loop repetitions
# reps = 1
methods = c("BH",                 # methods
            "Focused_BH_original",   
            "Focused_BH_permutation",
            "Yekutieli",
            "Structured_Holm")
q = 0.1                           # FDR control level
signal_strength_vals = seq(1,7,by=0.5)
# signal_strength_vals = 6
global_test = "Fisher"
B = 100       # number of repetitions
filter_name = "outer_nodes"
parameters = as_tibble(merge(1:reps, signal_strength_vals))
names(parameters) = c("rep", "signal_strength")
num_experiments = ceiling(reps*length(signal_strength_vals)/reps_per_experiment)
parameters$experiment = 1:num_experiments

### set input_mode (experiment, precomputation, num_experiments, num_precomputations)
if(exists("input_mode")){
  stopifnot(input_mode %in% c("experiment", "precomputation"))
  if(input_mode == "experiment"){
    stopifnot(exists("experiment_index"))
    parameters = parameters %>% dplyr::filter(experiment == experiment_index)
  }
  if(input_mode == "precomputation"){
    stopifnot(exists("b"))
  }
} else{
  args = commandArgs(trailingOnly = TRUE)
  num_args = length(args)
  stopifnot(num_args == 1)
  input_mode = args[1]
  stopifnot(input_mode %in% c("num_experiments", "num_precomputations"))
  if(input_mode == "num_experiments"){
    cat(num_experiments)
  }
  if(input_mode == "num_precomputations"){
    cat(B)
  }
}

### define the graph
if(input_mode %in% c("precomputation", "experiment")){
  # read in GO data
  graph_file = sprintf("%s/data/processed/ICD_graph.Rda", base_dir)
  load(graph_file)
  
  num_items = G$num_items
  m = G$m
  k_max = m
  
  # create adjacency matrix between genes and GO terms
  
  items_per_node = sapply(G$sets, length)
  i = integer(sum(items_per_node))
  j = integer(sum(items_per_node))
  index = 1
  for(node in 1:m){
    items_in_node = G$sets[[node]]
    i[index:(index+items_per_node[node]-1)] = items_in_node
    j[index:(index+items_per_node[node]-1)] = node
    index = index + items_per_node[node]
  }
  adj_matrix = sparseMatrix(i, j, dims = c(num_items, m))
  
  # Non-sparse equivalent code
  # adj_matrix = matrix(FALSE, num_items, m)
  # for(node in 1:m){
  #   items_in_node = G$sets[[node]]
  #   adj_matrix[items_in_node, node] = TRUE
  # }
  
  
  node_depths = get_depths(G$Pa)
  leaves = sapply(G$C, length) == 0
  num_deep_leaves = get_num_marked_descendants(G$C,  leaves & node_depths == 5)
  
  # sample 5 "anchor nodes" (ICD codes) with 10
  anchor_nodes = withSeed(sample(which(10 <= num_deep_leaves &
                                         num_deep_leaves <= 20 & node_depths == 3), 10), 1)

  # define the genes belonging to the anchor terms non-null
  nonnull_items = unique(unlist(lapply(G$sets[anchor_nodes], 
                                       function(set)(withSeed(sample(set, 5),1)))))
  
  # deduce which other GO terms are non-null
  nonnull_nodes = colSums(adj_matrix[nonnull_items,]) > 0
}