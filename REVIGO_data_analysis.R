rm(list = ls())
machine = "local"
source("setup.R")

# parameters
methods = c("Focused_BH", "Structured_Holm", "BH")
q = 0.1

cat(sprintf("Reading in GO data...\n"))
# read in GOrilla data
GO <- read_excel(sprintf("%s/data/raw/GO.xls", base_dir))

# read in GO data
reduced_graph_file = sprintf("%s/data/processed/reduced_graph.Rda", base_dir)
load(reduced_graph_file)
ids = names(G$gene_sets)
m = length(ids)
index = c()
for(i in 1:m){
  index[[ids[i]]] = i
}

P = rep(1, m)
significant_ids = match(GO$`GO Term`, ids)
cat(sprintf("%d out of %d of the significant GO terms are not in the older ontology.\n", 
            sum(is.na(significant_ids)), length(significant_ids)))
P[significant_ids[!is.na(significant_ids)]] = GO$`P-value`[!is.na(significant_ids)]
names(P) = ids

parents_indexed = sapply(G$Pa, function(parents_list)(unname(index[parents_list])))
children_indexed = sapply(G$C, function(children_list)(unname(index[children_list])))
sets = G$gene_sets

G = c()
G$m = m
G$Pa = parents_indexed
G$C = children_indexed

num_methods = length(methods)
rejections = matrix(FALSE, num_methods, m, dimnames = list(methods, ids))
rejections_filtered = matrix(FALSE, num_methods, m, dimnames = list(methods, ids))

# RUN STRUCTURED HOLM
cat(sprintf("Running Structured Holm...\n"))
G_SH = new("DAGstructure", parents = parents_indexed, children = children_indexed, sets = sets, twoway = FALSE)
if(min(P) >= q/m){
  R = rep(FALSE, m)
} else{
  R = structuredHolm(G_SH, alpha_max = q, pvalues=P)@rejected
}
R_F = run_filter(P, R, G, "REVIGO")
rejections["Structured_Holm",] = R
rejections_filtered["Structured_Holm",] = R_F

# RUN FOCUSED BH
cat(sprintf("Running Focused BH...\n"))
R = FocusedBH("REVIGO", P, G, q)
R_F = run_filter(P, R, G, "REVIGO")
rejections["Focused_BH",] = R
rejections_filtered["BH",] = R_F

# RUN BH
cat(sprintf("Running BH...\n"))
R = p.adjust(P, "fdr") <= q
R_F = run_filter(P, R, G, "REVIGO")
rejections["BH",] = R
rejections_filtered["BH",] = R_F

