suppressPackageStartupMessages(library(tidyverse))

######### SIMULATION PARAMETERS ###################
reps_per_experiment = 25               # number of repetitions per experiment
reps = 250                             # total number of repetitions
methods = c("BH",                      # methods to compare
            "Focused_BH_original",   
            "Focused_BH_permutation",
            "Structured_Holm")
q = 0.1                                # FDR control level
signal_strength_vals = seq(1,6,by=0.5) # signal strength
global_test = "Simes"                  # global test to aggregate item-level p-values
k_max = 350                            # maximum possible number of rejections
B = 100                                # number of permutations 
filter_name = "REVIGO"                 # filter

# put parameter combinations into a data frame
parameters = as_tibble(merge(1:reps, signal_strength_vals))
names(parameters) = c("rep", "signal_strength")
num_experiments = ceiling(reps*length(signal_strength_vals)/reps_per_experiment)
parameters$experiment = 1:num_experiments

######### SET INPUT MODE ###################
# Input mode is one of {"experiment", "precomputation", 
# "num_experiments", "num_precomputations"} and determines
# how this script behaves
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

######### DEFINE THE GRAPH AND NON-NULL STRUCTURE ###################
if(input_mode %in% c("precomputation", "experiment")){
  # read in GO graph
  GO_reduced_graph_100_file = sprintf("%s/data/processed/GO/GO_reduced_graph_100.Rda", base_dir)
  load(GO_reduced_graph_100_file)

  # create adjacency matrix between genes and GO terms
  num_items = G$num_items
  m = G$m
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