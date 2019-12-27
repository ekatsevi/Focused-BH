suppressPackageStartupMessages(library(tidyverse))

### set top-level parameters
reps_per_experiment = 25
reps = 100                         # Number of outer-loop repetitions
# reps = 1
methods = c("BH",                 # methods
            "Focused_BH_original",   
            "Focused_BH_permutation",
            "Yekutieli",
            # "Focused_BH_oracle",
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
  cat(sprintf("Reading in ICD graph...\n"))
  graph_file = sprintf("%s/data/processed/ICD_graph.Rda", base_dir)
  load(graph_file)
  
  gene_sets = G$gene_sets
  ids = names(G$C)
  m = G$m
  num_annotations = sapply(gene_sets, length)
  
  index = c()
  for(i in 1:m){
    index[[ids[i]]] = i
  }
  
  parents_indexed = sapply(G$Pa, function(parents_list)(unname(index[parents_list])))
  children_indexed = sapply(G$C, function(children_list)(unname(index[children_list])))
  
  G = c()
  G$m = m
  G$Pa = parents_indexed
  G$C = children_indexed
  
  if("Structured_Holm" %in% methods & input_mode == "experiment"){
    cat(sprintf("Creating DAG for Structured Holm...\n"))
    G_SH = new("DAGstructure", parents = parents_indexed, 
               children = children_indexed, sets = gene_sets, twoway = FALSE)
  }
  
  genes = unique(unlist(gene_sets))
  num_genes = length(genes)
  adj_matrix = matrix(FALSE, num_genes, m)
  rownames(adj_matrix) = genes
  colnames(adj_matrix) = ids
  for(id in ids){
    adj_matrix[gene_sets[[id]],id] = TRUE
  }
  
  anchor_terms = withSeed(sample(which(num_annotations == 5), 2), 1)
  nonnull_genes = unique(unlist(lapply(gene_sets[anchor_terms], function(set)(set))))
  nonnull_terms = colSums(adj_matrix[nonnull_genes,]) > 0
}