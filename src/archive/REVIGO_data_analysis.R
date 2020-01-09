# libraries
library(readxl)
library(ggplot2)
library(reshape2)
library(cherry)
library(tidyverse)

source("aux_DAG.R")

# parameters
dispensability_thresh = 0.5
q = 0.1
n_tested = 14398

# read in GOrilla data
GO <- read_excel("/home/ekatsevi/project-files/filtered_BH/DAG/data/raw/GO.xls")
n = nrow(GO)

# read in GO data
processed_data_file = "/home/ekatsevi/project-files/filtered_BH/DAG/data/processed/REVIGO_data.Rda"
load(processed_data_file)
m = length(ids)
index = c()
for(i in 1:m){
  index[[ids[i]]] = i
}

P = rep(1, m)
P[match(GO$`GO Term`, ids)] = GO$`P-value`

parents_indexed = sapply(parents, function(parents_list)(index[parents_list]))
children_indexed = sapply(children, function(children_list)(index[children_list]))
sets = vector("list", m)
names(sets) = ids

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

df_results = data.frame(ids[R], P[R])
SH_rejections_file = "SH_rejections.tsv"
write_tsv(x = df_results, path = SH_rejections_file, col_names = FALSE)

# rejected_indices = which(R)
# goterms = sprintf("%s %0.2e", ids[rejected_indices[1]], P[rejected_indices[1]])
# for(j in 2:length(rejected_indices)){
#   goterms = paste(goterms, sprintf("%s %0.2e", ids[rejected_indices[j]], P[rejected_indices[j]]), sep = "\n")
# }
# goterms = sprintf("\"\"\"\n%s\n\"\"\"", goterms)

REVIGO_path = "/home/ekatsevi/Dropbox/Research/Projects/HierTest/code/REVIGO"
command = sprintf("java -jar %s/RevigoStandalone.jar %s --stdout --cutoff=0.5", REVIGO_path, SH_rejections_file)

output = system(command, intern = TRUE)


# n_start = 135

FDP_hat = numeric(n_start)
FDP_hat_F = numeric(n_start)
R = numeric(n_start)
R_F = numeric(n_start)

n_start = 19
for(i in n_start:1){
  cat(sprintf("Processing first %d GO terms...\n", i))
  goterms = sprintf("%s %0.2e", GO[1,"GO Term"], GO[1,"P-value"])
  if(i > 1){
    for(j in 2:i){
      goterms = paste(goterms, sprintf("%s %0.2e", GO[j,"GO Term"], GO[j,"P-value"]), sep = "\n")
    }
  }
  goterms = sprintf("\"\"\"\n%s\n\"\"\"", goterms)
  
  command = sprintf("python scrape_revigo.py %s", goterms)
  output = system(command, intern = TRUE)
  
  num_terms = length(output)-2
  GO_terms = character(num_terms)
  dispensabilities = numeric(num_terms)
  for(j in 1:num_terms){
    GO_terms[j] = unlist(strsplit(output[j+1], split = ","))[1]
    dispensabilities[j] = as.numeric(unlist(strsplit(output[j+1], split = ","))[9])
  }
  discoveries = dispensabilities < dispensability_thresh
  n_discoveries = sum(discoveries) 
  
  R[i] = num_terms
  FDP_hat[i] = n_tested*GO[i,"P-value"]/num_terms
  R_F[i] = n_discoveries
  FDP_hat_F[i] = n_tested*GO[i,"P-value"]/n_discoveries
}

valid_idx = c(TRUE,R[2:n_start] > R[1:(n_start-1)])
n_start = 135
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


