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
reduced_graph_file = sprintf("%s/data/processed/GO_reduced_graph.Rda", base_dir)
load(reduced_graph_file)
ids = G$node_names
m = G$m

P = rep(1, m)
significant_ids = match(GO$`GO Term`, ids)
cat(sprintf("%d out of %d of the significant GO terms are not in the older ontology.\n", 
            sum(is.na(significant_ids)), length(significant_ids)))
P[significant_ids[!is.na(significant_ids)]] = GO$`P-value`[!is.na(significant_ids)]
names(P) = ids

num_methods = length(methods)
rejections = matrix(FALSE, num_methods, m, dimnames = list(methods, ids))
rejections_filtered = matrix(FALSE, num_methods, m, dimnames = list(methods, ids))

# RUN STRUCTURED HOLM
cat(sprintf("Running Structured Holm...\n"))
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
rejections_filtered["Focused_BH",] = R_F

# RUN BH
cat(sprintf("Running BH...\n"))
R = p.adjust(P, "fdr") <= q
R_F = run_filter(P, R, G, "REVIGO")
rejections["BH",] = R
rejections_filtered["BH",] = R_F

