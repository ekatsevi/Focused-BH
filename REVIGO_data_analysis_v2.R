source("setup.R")

# parameters
dispensability_thresh = 0.5
q = 0.1
n_tested = 14398

# read in GOrilla data
GO <- read_excel("/home/ekatsevi/project-files/filtered_BH/DAG/data/raw/GO.xls")
n = nrow(GO)

# read in GO data
reduced_graph_file = "/home/ekatsevi/project-files/filtered_BH/DAG/data/processed/reduced_graph.Rda"
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

parents_indexed = sapply(G$Pa, function(parents_list)(unname(index[parents_list])))
children_indexed = sapply(G$C, function(children_list)(unname(index[children_list])))
sets = G$gene_sets

G = c()
G$n = m
G$Pa = parents_indexed
G$C = children_indexed

# RUN STRUCTURED HOLM
q = 0.1
G_SH = new("DAGstructure", parents = parents_indexed, children = children_indexed, sets = sets, twoway = FALSE)

if(min(P) >= q/m){
  R = rep(FALSE, m)
} else{
  R = structuredHolm(G_SH, alpha_max = q, pvalues=P)@rejected
}
R = R & P < 0.5

df_results = data.frame(ids[R], P[R])
SH_rejections_file = "SH_rejections.tsv"
write_tsv(x = df_results, path = SH_rejections_file, col_names = FALSE)

rejected_indices = which(R)
goterms = sprintf("%s %0.2e", ids[rejected_indices[1]], P[rejected_indices[1]])
for(j in 2:length(rejected_indices)){
  goterms = paste(goterms, sprintf("%s %0.2e", ids[rejected_indices[j]], P[rejected_indices[j]]), sep = "\n")
}
goterms = sprintf("\"\"\"\n%s\n\"\"\"", goterms)

REVIGO_path = "/home/ekatsevi/Dropbox/Research/Projects/HierTest/code/REVIGO"
command = sprintf("java -jar %s/RevigoStandalone.jar %s --stdout --cutoff=%0.1f", 
                  REVIGO_path, SH_rejections_file, dispensability_thresh)
output = system(command, intern = TRUE)
matrix_output = do.call("rbind", lapply(output[3:length(output)], function(str)(unname(strsplit(str, split = "\t"))[[1]])))
df_output = as.data.frame(matrix_output[2:nrow(matrix_output),])
names(df_output) = matrix_output[1,]
rejections = as.character(df_output$term_ID[df_output$eliminated == 0])

ord = order(P)
n_start = length(significant_ids)

FDP_hat = numeric(n_start)
FDP_hat_F = numeric(n_start)
R = numeric(n_start)
R_F = numeric(n_start)

# n_start = 19
for(num_terms in n_start:1){
  cat(sprintf("Processing first %d GO terms...\n", num_terms))
  
  df_tmp = data.frame(ids[ord[1:num_terms]], P[ord[1:num_terms]])
  tmp_file = "rejections.tsv"
  write_tsv(x = df_tmp, path = tmp_file, col_names = FALSE)
  
  command = sprintf("java -jar %s/RevigoStandalone.jar %s --stdout --cutoff=%0.1f", 
                    REVIGO_path, tmp_file, dispensability_thresh)
  output = system(command, intern = TRUE)
  
  GO_terms = character(num_terms)
  dispensabilities = numeric(num_terms)
  eliminated = numeric(num_terms)
  for(j in 1:num_terms){
    GO_terms[j] = unlist(strsplit(output[j+3], split = "\t"))[2]
    dispensabilities[j] = as.numeric(unlist(strsplit(output[j+3], split = "\t"))[10])
    eliminated[j] = as.numeric(unlist(strsplit(output[j+3], split = "\t"))[12])
  }
  discoveries = dispensabilities < dispensability_thresh
  n_discoveries = sum(discoveries)
  
  R[num_terms] = num_terms
  FDP_hat[num_terms] = m*GO[num_terms,"P-value"]/num_terms
  R_F[num_terms] = n_discoveries
  FDP_hat_F[num_terms] = m*GO[num_terms,"P-value"]/n_discoveries
  if(FDP_hat_F[num_terms] <= q){
    break
  }
}

valid_idx = c(TRUE,R[2:n_start] > R[1:(n_start-1)])
# n_start = 135
n_valid = sum(valid_idx)

pvals = GO$"P-value"[1:n_start]


to_plot = data.frame(unlist(FDP_hat_F[valid_idx]), unlist(FDP_hat[valid_idx]), unlist(R_F[valid_idx]), unlist(R[valid_idx]))
names(to_plot) = c("FDP_hat_F", "FDP_hat", "R_F", "R")
to_plot$index = 1:n_valid
to_plot_reshaped = melt(to_plot, 5)

to_plot_F = to_plot[,c(1,3,5)]
names(to_plot_F) = c("FDP hat", "Number of discoveries", "Index")
to_plot_F$method = "Focused BH"

to_plot_U = to_plot[,c(2,4,5)]
names(to_plot_U) = c("FDP hat", "Number of discoveries", "Index")
to_plot_U$method = "BH"

to_plot_2 = rbind(to_plot_U, to_plot_F)

to_plot_3 = melt(to_plot_2, 3:4)

q = 0.1
value = q
variable = "FDP hat"
hline_data = data.frame(variable, value)

save(to_plot_3, file = "path_information.Rda")


thresh = max(which(FDP_hat[valid_idx] <= q))
thresh_F = max(which(FDP_hat_F[valid_idx] <= q))

p = ggplot(to_plot_3, aes(x = Index, y = value, group = method)) + facet_grid(variable ~ ., scales = "free") + 
  geom_point(aes(color = method)) + ylab("") + 
  geom_hline(data = hline_data, aes(yintercept = value), linetype = "dashed", color = "black") + 
  theme_bw(base_size = 14) + theme(legend.position = "bottom") + 
  scale_colour_manual(values = c("firebrick1", "dodgerblue"), name = "") + 
  geom_vline(xintercept = thresh, color = "firebrick1") + geom_vline(xintercept = thresh_F, color = "dodgerblue")
plot(p)
plot_filename = "/home/ekatsevi/Dropbox/Stanford/Research/Posters and Presentations/2018/Berkeley_Stats_Genomics_2018/figures/real_data.pdf"
ggsave(filename = plot_filename, plot = p, width = 4.5, height = 5, units = "in")


to_plot = data.frame(unlist(FDP_hat_F), unlist(FDP_hat), unlist(R_F), unlist(R))

plot(pvals[valid_idx], FDP_hat_F[valid_idx], pch = 19, col = "dodgerblue", log = "xy")
points(pvals[valid_idx], FDP_hat[valid_idx], pch = 19, col = "firebrick1")
abline(h = 0.1, lty = 2)

plot(1:n_valid, R[valid_idx], pch = 19, col = "firebrick1")
points(1:n_valid, R_F[valid_idx], pch = 19, col = "dodgerblue")


