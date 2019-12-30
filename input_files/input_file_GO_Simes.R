suppressPackageStartupMessages(library(tidyverse))

### set top-level parameters
reps_per_experiment = 25
reps = 250                         # Number of outer-loop repetitions
# reps = 1
methods = c("BH",                 # methods
            "Focused_BH_original",   
            "Focused_BH_permutation",
            "Structured_Holm")
q = 0.1                           # FDR control level
signal_strength_vals = seq(1,6,by=0.5)
# signal_strength_vals = 6
global_test = "Simes"
k_max = 350   # maximum p-value threshold considered by permutation approach
B = 100       # number of repetitions
filter_name = "REVIGO"
parameters = as_tibble(merge(1:reps, signal_strength_vals))
names(parameters) = c("rep", "signal_strength")
num_experiments = ceiling(reps*length(signal_strength_vals)/reps_per_experiment)
parameters$experiment = 1:num_experiments

### set input_mode (experiment, precomputation, num_experiments, num_precomputations)
# TBD: write a function that does the following
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
  GO_reduced_graph_100_file = sprintf("%s/data/processed/GO_reduced_graph_100.Rda", base_dir)
  load(GO_reduced_graph_100_file)

  num_items = G$num_items
  m = G$m
  # create adjacency matrix between genes and GO terms
  adj_matrix = matrix(FALSE, num_items, m)
  for(node in 1:m){
    genes_in_node = G$sets[[node]]
    adj_matrix[genes_in_node, node] = TRUE
  }
  items_per_node = sapply(G$sets, length)
  
  # sample 2 "anchor nodes" (GO terms) with 5 annotations each
  anchor_nodes = withSeed(sample(which(items_per_node == 5), 2), 1)
  
  # define the genes belonging to the anchor terms non-null
  nonnull_items = unique(unlist(G$sets[anchor_nodes]))
  
  # deduce which other GO terms are non-null
  nonnull_nodes = colSums(adj_matrix[nonnull_items,]) > 0
}