library(tidyverse)

### set top-level parameters
reps = 10                         # Number of outer-loop repetitions
methods = c("BH",                 # methods
            "Focused_BH_original",   
            "Focused_BH_permutation",
            "Structured_Holm")
q = 0.1                           # FDR control level
signal_strength_vals = c(1,2,3)
global_test = "Fisher"
t_max = 0.02  # maximum p-value threshold considered by permutation approach
B = 100       # number of repetitions
filter_name = "REVIGO"
parameters = tibble(signal_strength_vals)
names(parameters) = "signal_strength"

### set mode (experiment, precomputation, num_experiments, num_precomputations)
args = commandArgs(trailingOnly = TRUE)
num_args = length(args)
stopifnot(num_args <= 1)
if(num_args == 0){
  stopifnot(exists("mode"))
  stopifnot("mode" %in% c("experiment", "precomputation"))
  if(mode == "experiment"){
    stopifnot(exists("experiment_index"))
  }
  if(mode == "precomputation"){
    stopifnot(exists("b"))
  }
}
if(num_args == 1){
  mode = args[1]
  stopifnot("mode" %in% c("num_experiments", "num_precomputations"))
  if(mode == "num_experiments"){
    cat(nrow(parameters))
  }
  if(mode == "num_precomputations"){
    cat(B)
  }
}

### define the graph
if(mode %in% c("precomputation", "experiment")){
  # read in GO data
  cat(sprintf("Reading in GO data...\n"))
  # reduced_graph_50_file = sprintf("%s/data/processed/reduced_graph_50.Rda", base_dir)
  # load(reduced_graph_50_file)
  reduced_graph_file = sprintf("%s/data/processed/reduced_graph.Rda", base_dir)
  load(reduced_graph_file)
  
  ids = names(G$gene_sets)
  m = length(ids)
  index = c()
  for(i in 1:m){
    index[[ids[i]]] = i
  }
  
  parents_indexed = sapply(G$Pa, function(parents_list)(unname(index[parents_list])))
  children_indexed = sapply(G$C, function(children_list)(unname(index[children_list])))
  sets = G$gene_sets
  
  G = c()
  G$m = m
  G$Pa = parents_indexed
  G$C = children_indexed
  
  if("Structured_Holm" %in% methods){
    cat(sprintf("Creating DAG for Structured Holm...\n"))
    G_SH = new("DAGstructure", parents = parents_indexed, 
               children = children_indexed, sets = sets, twoway = FALSE)
  }
  
  genes = unique(unlist(sets))
  num_genes = length(genes)
  adj_matrix = matrix(FALSE, num_genes, m)
  rownames(adj_matrix) = genes
  colnames(adj_matrix) = ids
  for(id in ids){
    adj_matrix[sets[[id]],id] = TRUE
  }
  
  anchor_terms = withSeed(sample(which(num_annotations == 5), 4), 1)
  nonnull_genes = unique(unlist(lapply(sets[anchor_terms], function(set)(set))))
  nonnull_terms = colSums(adj_matrix[nonnull_genes,]) > 0  
  
  # mean(nonnull_terms)
  # P = numeric(m)
  # P[!nonnull_terms] = 1
  # R_F = run_filter(P, nonnull_terms, G, "REVIGO")
  # hist(colSums(adj_matrix[nonnull_genes,nonnull_terms]), breaks = 10)
  # hist(colSums(adj_matrix[nonnull_genes,nonnull_terms])/colSums(adj_matrix[,nonnull_terms]), breaks = 20)
}