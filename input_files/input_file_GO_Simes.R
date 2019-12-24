suppressPackageStartupMessages(library(tidyverse))

### set top-level parameters
reps = 2                         # Number of outer-loop repetitions
methods = c("BH",                 # methods
            "Focused_BH_original",   
            "Focused_BH_permutation",
            "Focused_BH_oracle",
            "Structured_Holm")
q = 0.1                           # FDR control level
signal_strength_vals = c(1,2,3,4,5)
global_test = "Simes"
k_max = 350  # maximum p-value threshold considered by permutation approach
B = 100       # number of repetitions
filter_name = "REVIGO"
parameters = tibble(signal_strength_vals)
names(parameters) = "signal_strength"

### set input_mode (experiment, precomputation, num_experiments, num_precomputations)
if(exists("input_mode")){
  stopifnot(input_mode %in% c("experiment", "precomputation"))
  if(input_mode == "experiment"){
    stopifnot(exists("experiment_index"))
    signal_strength = parameters$signal_strength[experiment_index]
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
    cat(nrow(parameters))
  }
  if(input_mode == "num_precomputations"){
    cat(B)
  }
}

### define the graph
if(input_mode %in% c("precomputation", "experiment")){
  # read in GO data
  cat(sprintf("Reading in GO data...\n"))
  # reduced_graph_50_file = sprintf("%s/data/processed/reduced_graph_50.Rda", base_dir)
  # load(reduced_graph_50_file)
  reduced_graph_file = sprintf("%s/data/processed/reduced_graph.Rda", base_dir)
  load(reduced_graph_file)
  
  gene_sets = G$gene_sets
  to_remove = names(which(sapply(G$gene_sets, length) > 100))
  G = subset_graph(to_remove, G)
  ids = names(G$C)
  m = length(ids)
  gene_sets = gene_sets[ids]
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
  
  anchor_terms = withSeed(sample(which(num_annotations == 5), 3), 1)
  nonnull_genes = unique(unlist(lapply(gene_sets[anchor_terms], function(set)(set))))
  nonnull_terms = colSums(adj_matrix[nonnull_genes,]) > 0
  
  # sum(nonnull_terms)
  # mean(nonnull_terms)
  # P = numeric(m)
  # P[!nonnull_terms | num_annotations > 100] = 1
  # R_F = run_filter(P, nonnull_terms & num_annotations <= 100, G, "REVIGO")
  # hist(colSums(adj_matrix[nonnull_genes,nonnull_terms]), breaks = 10)
  # hist(colSums(adj_matrix[nonnull_genes,nonnull_terms])/colSums(adj_matrix[,nonnull_terms]), breaks = 20)
}