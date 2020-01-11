suppressPackageStartupMessages(library(tidyverse))

######### SIMULATION PARAMETERS ###################
reps_per_experiment = 250              # repetitions per experiment
reps = 250                             # total number of repetitions
methods = c("BH",                      # methods to compare
            "Focused_BH_original",   
            "Yekutieli",
            "Structured_Holm")
q = 0.1                                # FDR control level
signal_strength_vals = seq(1,6,by=0.5) # signal strength
global_test = "Fisher"                 # global test to aggregate item-level p-values
B = 100                                # number of permutations
filter_name = "outer_nodes"            # filter

# put parameter combinations into a data frame
parameters = as_tibble(merge(1:reps, signal_strength_vals))
names(parameters) = c("rep", "signal_strength")
num_experiments = ceiling(reps*length(signal_strength_vals)/reps_per_experiment)
parameters$experiment = 1:num_experiments

######### SET INPUT MODE ###################
# Input mode is one of {"experiment", "precomputation", 
# "num_experiments", "num_precomputations"} and determines
# how this script behaves
if(!exists("input_mode")){
  args = commandArgs(trailingOnly = TRUE)
  num_args = length(args)
  stopifnot(num_args == 1)
  input_mode = args[1]
}
if(input_mode == "experiment"){
  stopifnot(exists("experiment_index"))
  parameters = parameters %>% dplyr::filter(experiment == experiment_index)
}
if(input_mode == "precomputation"){
  stopifnot(exists("b"))
}
if(input_mode == "num_experiments"){
  cat(num_experiments)
}
if(input_mode == "num_precomputations"){
  cat(B)
}

######### DEFINE THE GRAPH AND NON-NULL STRUCTURE ###################
if(input_mode %in% c("precomputation", "experiment")){
  # read in GO data
  graph_file = sprintf("%s/data/processed/biobank/ICD_graph.Rda", base_dir)
  load(graph_file)
  
  
  # create adjacency matrix between genes and GO terms
  num_items = G$num_items
  m = G$m
  k_max = m
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
  
  # sample 10 "anchor nodes" (ICD codes) with 5 items each, in different subtrees
  node_depths = get_depths(G$Pa)
  leaves = sapply(G$C, length) == 0
  num_leaves = get_num_marked_descendants(G$C,  leaves)
  root_nodes = get_root_nodes(G$Pa)
  roots = which(node_depths == 1)
  nonnull_roots = withSeed(sample(roots, 10),1)
  anchor_nodes = unlist(sapply(nonnull_roots, 
                               function(nonnull_root)(withSeed(sample(which(root_nodes == nonnull_root & num_leaves == 5),1),1))))
  
  # define the genes belonging to the anchor terms non-null
  nonnull_items = unlist(G$sets[anchor_nodes])
  
  # deduce which other GO terms are non-null
  nonnull_nodes = colSums(adj_matrix[nonnull_items,]) > 0
}