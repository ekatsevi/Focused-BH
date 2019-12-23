rm(list = ls())

#allele = "B_2705"    # which HLA allele we are looking at
allele = "DRB1_301"
q = 0.1           # target FDR level

setwd("~/Dropbox/Stanford/Research/Projects/HierTest/code/code_DAG/src/UKBiobank_PheWAS")
source("aux_tree.R")
library(readr)
library(Matrix)
options(tibble.print_max = Inf) 

######### RAW DATA FILENAMES ###################
tree_file = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank/coding6.tsv"
phen_file = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank/self_reported_disease.tab"
geno_file = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank/HLA.tab"
allele_file = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank/haplotype_names.txt"

######### PROCESS TREE DATA #############

# read tree
tree = read_tsv(tree_file)
# remove last row (coding 99999, unclassifiable)
tree = tree[1:(nrow(tree)-1),]
m = nrow(tree)
id2idx = list()
for(j in 1:m){
  id2idx[as.character(tree$node_id[j])] = j
}
C = vector("list", m)
Pa = vector("list", m)
for(j in 1:m){
  node_id = tree$node_id[j]
  parent_id = tree$parent_id[j]
  if(parent_id != 0){
    parent_idx = id2idx[[as.character(parent_id)]]
    Pa[[j]] = c(Pa[[j]], parent_idx)
    C[[parent_idx]] = c(C[[parent_idx]], j)
  }
}

G = c()
G$m = m
G$C = C
G$Pa = Pa
G$depths = get_depths(Pa)


######### PROCESS PHENOTYPE DATA #############
SR_affected_file = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank/SR_affected.Rda"
if(file.exists(SR_affected_file)){
  load(SR_affected_file)
} else{
  SR = read_tsv(phen_file, col_types = cols(.default = col_integer()))
  n = nrow(SR)
  p = ncol(SR)
  diseases = unique(unlist(sapply(1:p, function(j)(unique(SR[[j]])))))
  diseases = diseases[!is.na(diseases)]
  m = length(diseases)
  affected = vector("list", m)
  names(affected) = diseases
  for(i in 1:n){
    if(i %% 10000 == 0){
      print(i)
    }
    for(j in 1:p){
      disease = SR[i,j]
      if(!is.na(disease)){
        affected[[as.character(disease)]] = c(affected[[as.character(disease)]], i)
      } else{
        break
      }
    }
  }
  new_names = sapply(names(affected), function(name)(tree$node_id[tree$coding == as.numeric(name)]))
  names(affected) = new_names
  for(j in which(tree$selectable == "N")){
    print(j)
    affecteds = c()
    to_explore = G$C[[j]] # should be non-empty
    while(length(to_explore) > 0){
      print(to_explore)
      k = to_explore[1]
      if(tree$selectable[k] == "Y"){
        affecteds = c(affecteds, affected[[as.character(tree$node_id[k])]]) 
      } else{
        to_explore = c(to_explore, G$C[[k]])
      }
      if(length(to_explore) == 1){
        to_explore = numeric(0)
      } else{
        to_explore = to_explore[2:length(to_explore)]
      }
    }
    affected[[as.character(tree$node_id[j])]] = affecteds
  }
  save(affected, file = SR_affected_file)
  # somewhere here, remove the 99999 field ("unclassifiable")
}

######### PROCESS HLA DATA #############
HLA_split_file = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank/HLA_split.tsv"

if(file.exists(HLA_split_file)){
  HLA_split_df = read_tsv(HLA_split_file)
} else{
  HLA = read_tsv(geno_file)
  HLA_names = read_tsv(allele_file, col_names = FALSE)
  n = nrow(HLA)
  p = length(strsplit(HLA$f.22182.0.0[1], split = ",")[[1]])
  HLA_split = matrix(0, n, p)
  for(i in 1:n){
    if(i %% 10000 == 0){
      print(i)
    }
    HLA_split[i,] = as.numeric(strsplit(HLA$f.22182.0.0[i], split = ",")[[1]])
  }
  HLA_split_df = data.frame(HLA_split)
  names(HLA_split_df) = HLA_names$X1
  write_tsv(HLA_split_df, path = HLA_split_file, col_names = TRUE)
}


######### CREATE P-VALUES ###############

genotypes = round(HLA_split_df[[allele]])
P = numeric(m)
for(k in 1:m){
  if(k %% 10 == 0){
    print(k)
  }
  affecteds = affected[[as.character(tree$node_id[k])]]
  aff = c(sum(genotypes[affecteds] == 0, na.rm = TRUE), 
          sum(genotypes[affecteds] == 1, na.rm = TRUE), 
          sum(genotypes[affecteds] == 2, na.rm = TRUE))
  unaff = c(sum(genotypes == 0, na.rm = TRUE), 
            sum(genotypes == 1, na.rm = TRUE), 
            sum(genotypes == 2, na.rm = TRUE)) - aff
  cont_table = data.frame(aff, unaff)
  P[k] = chisq.test(cont_table)$p.value
}
P[is.na(P)] = 1

####### CARRY OUT THE TESTS ###############

# BH

BH_rejections = p.adjust(P, "fdr") <= q
sum(BH_rejections)
sum(get_outer_nodes(BH_rejections, G))

ancestors = get_ancestors(BH_rejections, G)
G_anc = subset_graph(ancestors,G)

#not_selectable = tree$selectable == "N"

vertex_shapes = "circle"
frame_colors = NULL
fill_colors = rep("white", G_anc$m)
fill_colors[BH_rejections[ancestors]] = "red"
#fill_colors[not_selectable[ancestors]] = "gray"
#vertex_labels = which(ancestors)
vertex_labels = NULL
vertex_size = 5
plot_tree(G_anc, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)



# Focused BH

gamma = function(p)(m*p)
num_rejections = get_num_outer_nodes(P, G)
P_aug = sort(c(0,P), decreasing = TRUE)
FDP_hat = gamma(P_aug)/num_rejections
FDP_hat[is.na(FDP_hat)] = 0
FDP_hat[G$m+1] = 0
idx = min(which(FDP_hat <= q))
#candidates = which(FDP_hat <= q)
#idx = candidates[which.max(num_rejections[candidates])]
R = P <= P_aug[idx]
sum(R)
sum(get_outer_nodes(R, G))

