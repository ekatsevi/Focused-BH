set.seed(1234)                    # for replicability

reps = 10                         # Number of outer-loop repetitions
methods = c("BH", "Focused_BH",   # methods
            "Structured_Holm")
q = 0.1                           # FDR control level
signal_strength_vals = c(1,2,3)

nonnull_genes = 1:50

k_max = 35 # maximum number of rejections considered
B = 100      # number of repetitions

filter_name = "REVIGO"
# base_log_odds = qlogis(0.05)

# define the graph
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

cat(sprintf("Creating DAG for Structured Holm...\n"))
G_SH = new("DAGstructure", parents = parents_indexed, children = children_indexed, sets = sets, twoway = FALSE)

genes = unique(unlist(sets))
num_genes = length(genes)
adj_matrix = matrix(FALSE, num_genes, m)
rownames(adj_matrix) = genes
colnames(adj_matrix) = ids
for(id in ids){
  adj_matrix[sets[[id]],id] = TRUE
}

nonnull_terms = colSums(adj_matrix[nonnull_genes,]) > 0  

parameters = tibble(signal_strength_vals)
names(parameters) = "signal_strength"

if(exists("experiment_index")){
  # set parameters
  signal_strength = parameters$signal_strength[experiment_index]
} else{
  # return number of parameters (for submission script)
  cat(nrow(parameters))  
}