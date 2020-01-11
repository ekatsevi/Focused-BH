######################################################
# Demonstrate the application of FocusedBH()
#
# In this file, we conduct a simulated phenome-wide 
# association study with the tree structure 
######################################################

# set up workspace
rm(list = ls())
set.seed(1234) # for replicability
source("FocusedBH.R")
source("aux.R")
library(readr)
library(Matrix)

# analysis parameters
q = 0.1                     # target FDR level
tree_file = "coding6.tsv"   # file with tree structure downloaded from UK Biobank 
                            # (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)
allele_frequency = 0.1      # simulated allele frequency for SNP
disease_frequency = 0.01    # baseline disease frequency 
                            # (the same for all leaf node diseases, for simplicity)
effect_size = 0.5           # effect size of non-reference allele at causal SNPs
n = 10000                   # sample size
num_associated_leaves = 10  # number of leaf nodes associated with SNP

# read tree structure from file and create graph object G
tree = read_tsv(tree_file)
m = nrow(tree) # number of nodes in the tree
id2idx = list()
for(j in 1:m){
  id2idx[as.character(tree$node_id[j])] = j
}
C = vector("list", m)  # the jth entry is the list of children of node j
Pa = vector("list", m) # the jth entry is the list of parents of node j
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

# simulate genotype data based on allele frequency
genotypes = rbinom(n = n, size = 2, prob = allele_frequency)

# simulate phenotype data
leaves = sapply(G$C, length) == 0                  # set of all leaves
associated_leaves = logical(m)                     # set of leaves associated with SNP
associated_leaves[which(leaves)[1:num_associated_leaves]] = TRUE
unassociated_leaves = leaves & !associated_leaves  # set of leaves not associated with SNP

affected = matrix(0, n, m) # affected[i,j] = 1 if individual i is affected by disease j; 
                           #               = 0 otherwise
# sample diseases from unassociated leaves from baseline disease frequency
for(leaf in which(unassociated_leaves)){ 
  affected[,leaf] = rbinom(n = n, size = 1, prob = disease_frequency)  
}
# sample diseases from associated leaves based on logistic model 
for(leaf in which(associated_leaves)){
  baseline = log(disease_frequency/(1-disease_frequency))
  affected[,leaf] = rbinom(n = n, size = 1, prob = 1/(1 + exp(-(baseline+ genotypes*effect_size))))
}
# individuals affected by a non-leaf disease are those that are 
# affected by any of the descendant leaves
for(node in which(!leaves)){
  node_logical = logical(m)
  node_logical[node] = TRUE
  leaf_descendants = get_descendants(node_logical, G) & leaves
  affected[,node] = rowSums(affected[,leaf_descendants,drop = F]) > 0
}

# calculate p-values by correlating genotype with phenotype
P = sapply(1:m, function(j)(cor.test(affected[,j], genotypes)$p.value))

# apply Focused BH
R = FocusedBH(P, G, q)
outer_nodes = get_outer_nodes(R, G)

# print results
cat(sprintf("Focused BH rejected %d nodes, including %d outer nodes.\n", 
            sum(R), sum(outer_nodes)))
node = which(R)
name = tree$meaning[R]
outer = outer_nodes[node]
results = data.frame(node, name, outer)
print(results)

# plot results
plot_rejections(R, G)
