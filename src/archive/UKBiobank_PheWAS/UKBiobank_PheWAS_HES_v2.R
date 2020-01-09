rm(list = ls())

allele = "B_2705"    # which HLA allele we are looking at
#allele = "DRB1_301"
q = 0.05
min_affected = 50
data_dir = "/home/ekatsevi/Local/large_file_storage/filtered_BH/biobank"

setwd("/home/ekatsevi/Dropbox/Research/Projects/2018/HierTest/code/UKBiobank_PheWAS")
source("aux_tree.R")
library(readr)
library(Matrix)
library(reshape2)
library(ggplot2)
options(tibble.print_max = Inf) 

######### RAW DATA FILENAMES ###################
tree_file = sprintf("%s/coding19.tsv", data_dir)
phen_file = sprintf("%s/diagnosed_disease.tab", data_dir)
geno_file = sprintf("%s/HLA.tab", data_dir)
allele_file = sprintf("%s/haplotype_names.txt", data_dir)

######### PROCESS TREE DATA #############

# read tree
tree = read_tsv(tree_file)
m = nrow(tree)

# create graph structure
graph_file = sprintf("%s/graph.Rda", data_dir)
if(!file.exists(graph_file)){
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
  
  save(G, file = graph_file)
}
load(graph_file)


######### PROCESS PHENOTYPE DATA #############
HES_affected_file = sprintf("%s/HES_affected.Rda", data_dir)

if(!file.exists(HES_affected_file)){
  HES = read_tsv(phen_file, col_types = cols(.default = col_character()))
  n = nrow(HES)
  p = ncol(HES)
  diseases = unique(unlist(sapply(1:p, function(j)(unique(HES[[j]])))))
  diseases = diseases[!is.na(diseases)]
  m0 = length(diseases)
  affected = vector("list", m0)
  names(affected) = diseases
  for(i in 1:n){
    if(i %% 10000 == 0){
      print(i)
    }
    for(j in 1:p){
      disease = HES[i,j]
      if(!is.na(disease)){
        affected[[as.character(disease)]] = c(affected[[as.character(disease)]], i)
      } else{
        break
      }
    }
  }
  new_names = sapply(names(affected), function(name)(tree$node_id[tree$coding == name]))
  names(affected) = new_names
  for(j in which(tree$selectable == "N")){
    affecteds = c()
    to_explore = G$C[[j]] # should be non-empty
    while(length(to_explore) > 0){
      #      print(to_explore)
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
  save(affected, file = HES_affected_file)
}

load(HES_affected_file)
num_affected = sapply(1:m, function(k)(length(affected[[as.character(tree$node_id[k])]])))

######### PROCESS HLA DATA #############
HLA_split_file = sprintf("%s/HLA_split.tsv", data_dir)

if(!file.exists(HLA_split_file)){
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

HLA_split_df = read_tsv(HLA_split_file)

######### CREATE P-VALUES ###############

pvals_filename = sprintf("%s/pvals_HES_%s_m%d.tsv", data_dir, allele, min_affected)
if(!file.exists(pvals_filename)){
  genotypes = round(HLA_split_df[[allele]])
  P = numeric(m)
  P[num_affected < min_affected] = 1
  for(k in which(num_affected >= min_affected)){
    # if(k %% 100 == 0){
    #   print(k)
    # }
    print(k)
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
  write_tsv(data.frame(P), path = pvals_filename)
}

P = read_tsv(pvals_filename)
P = P$P
stopifnot(!any(is.na(P)))
#P[is.na(P)] = 1

####### CARRY OUT THE TESTS ###############

# BH

P[num_affected < min_affected] = 1
m_tested = sum(num_affected >= min_affected)

BH_rejections = p.adjust(P, "fdr") <= q*m/m_tested
BH_outer_nodes = get_outer_nodes(BH_rejections, G)
sum(BH_rejections)
sum(BH_outer_nodes)

ancestors = get_ancestors(BH_rejections, G)
G_anc = subset_graph(ancestors,G)

vertex_shapes = "circle"
frame_colors = NULL
fill_colors = rep("white", G_anc$m)
fill_colors[BH_rejections[ancestors]] = "gray35"
fill_colors[BH_outer_nodes[ancestors]] = "red"
vertex_labels = which(ancestors)
#vertex_labels = NULL
vertex_size = 5
plot_tree(G_anc, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)

# Focused BH
gamma = function(p)(m_tested*p)
num_rejections = get_num_outer_nodes(P, G)
P_aug = sort(c(0,P), decreasing = TRUE)
FDP_hat = gamma(P_aug)/num_rejections
FDP_hat[is.na(FDP_hat)] = 0
FDP_hat[G$m+1] = 0
idx = min(which(FDP_hat <= q))
#candidates = which(FDP_hat <= q)
#idx = candidates[which.max(num_rejections[candidates])]
FBH_rejections = P <= P_aug[idx]
FBH_outer_nodes = get_outer_nodes(FBH_rejections, G)
sum(FBH_rejections)
sum(FBH_outer_nodes)

# Yekutieli
library(structSSI)   # to perform Yekutieli test
library(igraph)      # to perform Yekutieli test
library(ape)         # to perform Yekutieli test

G_to_test = subset_graph(num_affected >= min_affected,G)
roots = which(sapply(G_to_test$Pa, length) == 0)
m_tested = G_to_test$m
G_to_test$m = m_tested + 1
G_to_test$C[[m_tested+1]] = roots
G_to_test$Pa[[m_tested+1]] = integer(0)
for(root in roots){
  G_to_test$Pa[[root]] = m_tested + 1
}
#G_to_test$depths = c(G_to_test$depths + 1, 1)
m_tested = m_tested + 1

G_Yekutieli = t(do.call(cbind, 
                        lapply(1:G_to_test$m, 
                               function(i)(if(length(G_to_test$C[[i]]) > 0) sapply(G_to_test$C[[i]], 
                                                                           function(j)(c(as.character(i), as.character(j)))) else matrix(0,2,0)))))
q_Yekutieli = q/(2*max(G_to_test$depths))
unadj.p.values = c(P[num_affected >= min_affected], 0)
names(unadj.p.values) = sapply(1:m_tested, function(i)(as.character(i)))
hyp.tree = hFDR.adjust(unadj.p.values, G_Yekutieli, q_Yekutieli)
Yekutieli_rejections = hyp.tree@p.vals[[2]] <= q_Yekutieli
Yekutieli_rejections[is.na(Yekutieli_rejections)] = FALSE
Yekutieli_rejections = Yekutieli_rejections[1:(m_tested-1)]
Y_rejections = logical(m)
Y_rejections[num_affected >= min_affected] = Yekutieli_rejections
sum(Y_rejections)
Y_outer_nodes = get_outer_nodes(Y_rejections, G)


# Structured Holm
library(cherry)      # to run the Structured Holm procedure (Meijer and Goeman, 2015b)
G_to_test = subset_graph(num_affected >= min_affected,G)
leaves = sapply(G_to_test$C, length) == 0
groups = vector("list", G_to_test$m)
for(j in 1:G_to_test$m){
  R_delta = logical(G_to_test$m)
  R_delta[j] = TRUE
  groups[[j]] = which(get_descendants(R_delta, G_to_test) & leaves)
}
G_to_test$groups = groups

get_descendants = function(R, G){
  m = G$m
  C = G$C
  traverse = function(node, descendant){
    if(!descendant[node]){
      descendant[node] = TRUE
      for(child in C[[node]]){
        descendant = traverse(child, descendant)
      }
    }
    return(descendant)
  }
  descendant = logical(m)
  for(node in which(R)){
    descendant = traverse(node, descendant)
  }
  return(descendant)
}

construct_SH = function (G) 
{
  graph = construct(G$groups)
  
  if(length(graph@parents) != G$m){
    graph = new("DAGstructure", parents = G$Pa, children = G$C, sets = G$groups, twoway = FALSE)
  }
  
  return(graph)
}
G_Structured_Holm = construct_SH(G_to_test)
SH_rejections_tmp = structuredHolm(G_Structured_Holm, alpha_max = q, 
                   pvalues=P[num_affected >= min_affected])@rejected
SH_rejections = logical(m)
SH_rejections[num_affected >= min_affected] = SH_rejections_tmp
SH_outer_nodes = get_outer_nodes(SH_rejections, G)
sum(SH_rejections)
sum(SH_outer_nodes)

######################################################################

# Plot progress of each algorithm
tested = num_affected >= min_affected
R = 1:m_tested
R_F = num_rejections[1:nrow(tree)][sort(P, decreasing = TRUE) < 1][m_tested:1]
FDP_hat = m_tested*sort(P[tested])/R
FDP_hat_F = m_tested*sort(P[tested])/R_F

to_plot = data.frame(FDP_hat_F, FDP_hat, R_F, R)
names(to_plot) = c("FDP_hat_F", "FDP_hat", "R_F", "R")
to_plot$index = 1:m_tested
to_plot_reshaped = melt(to_plot, 5)

to_plot_F = to_plot[,c(1,3,5)]
names(to_plot_F) = c("FDP hat", "Number of discoveries", "Index")
to_plot_F$method = "Focused BH"

to_plot_U = to_plot[,c(2,4,5)]
names(to_plot_U) = c("FDP hat", "Number of discoveries", "Index")
to_plot_U$method = "BH"

to_plot_2 = rbind(to_plot_U, to_plot_F)

to_plot_3 = melt(to_plot_2, 3:4)
to_plot_3$variable = factor(to_plot_3$variable, levels = c("Number of discoveries", "FDP hat"))

value = q
variable = "FDP hat"
hline_data = data.frame(variable, value)

#save(to_plot_3, file = "path_information.Rda")


thresh = max(which(FDP_hat <= q))
thresh_F = max(which(FDP_hat_F <= q))

p = ggplot(to_plot_3[to_plot_3$Index <= 75,], aes(x = Index, y = value, group = method)) + facet_grid(variable ~ ., scales = "free") + 
  geom_line(aes(color = method), size = 1.5) + ylab("") + 
  geom_hline(data = hline_data, aes(yintercept = value), linetype = "dashed", color = "black") + 
  theme_bw(base_size = 14) + theme(legend.position = "bottom") + 
  scale_colour_manual(values = c("firebrick1", "dodgerblue"), name = "") + 
  geom_vline(xintercept = thresh, color = "firebrick1") + geom_vline(xintercept = thresh_F, color = "dodgerblue")
plot(p)
plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/FBH_progress.pdf"
ggsave(filename = plot_filename, plot = p, width = 9, height = 6, units = "in")



######################################################
ord = order(P[BH_rejections])


vertex_shapes = "circle"
frame_colors = rep("black", G_anc$m)
frame_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "red"
fill_colors = rep("white", G_anc$m)
fill_colors[BH_rejections[ancestors]] = "gray45"
fill_colors[BH_outer_nodes[ancestors]] = "red"
fill_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "white"
vertex_labels = rep("", G_anc$m)
vertex_labels[BH_rejections[ancestors]] = invPerm(ord)
#vertex_labels = NULL
vertex_size = 5
plot_tree(G_anc, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)

final_results = cbind(tree[BH_rejections,2], P[BH_rejections])
colnames(final_results) = c("ICD-10 term", "p-value")
final_results$outer_node = FALSE
final_results$outer_node[BH_outer_nodes[BH_rejections]] = TRUE
final_results$FBH= FALSE
final_results$FBH[FBH_rejections[BH_rejections]] = TRUE
final_results = final_results[ord,]
rownames(final_results) = NULL
write_tsv(final_results, path = "B2705_rejections.txt")

##############################################################################

# Bird's eye view
vertex_shapes = "square"
frame_colors = rep("black", G_anc$m)
frame_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "darkcyan"
fill_colors = rep("white", G_anc$m)
fill_colors[BH_rejections[ancestors]] = "gray45"
fill_colors[BH_outer_nodes[ancestors]] = "darkcyan"
fill_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "white"
#vertex_labels = rep("", G_anc$m)
vertex_labels = rep(1:G_anc$m, G_anc$m)
#vertex_labels[BH_rejections[ancestors]] = invPerm(ord)
#vertex_labels = NULL
vertex_size = 5

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/output_ICD_tree.pdf"
pdf(file = plot_filename, width = 16, height = 10)
par(mar=c(0,0,0,0)) # make smaller margins
plot_tree(G_anc, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)
dev.off()
command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)

##############################################################################

# Bird's eye view, focus AS
R = logical(G$m)
R[1756] = TRUE
descendants = get_descendants(R, G)
to_plot = descendants & ancestors

vertex_shapes = "square"
frame_colors = rep("black", G_anc$m)
frame_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "darkcyan"
frame_colors[!to_plot[ancestors]] = "lightgray"
fill_colors = rep("white", G_anc$m)
fill_colors[BH_rejections[ancestors]] = "gray45"
fill_colors[BH_outer_nodes[ancestors]] = "darkcyan"
fill_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "white"
fill_colors[!to_plot[ancestors]] = "lightgray"
vertex_labels = rep("", G_anc$m)
#vertex_labels[BH_rejections[ancestors]] = invPerm(ord)
#vertex_labels = NULL
vertex_size = 5

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/output_ICD_tree_focus_AS.pdf"
pdf(file = plot_filename, width = 16, height = 10)
par(mar=c(0,0,0,0)) # make smaller margins
plot_tree(G_anc, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)
dev.off()
command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)
#title("Diseases of the musculoskeletal system and connective tissue")

##############################################################################

# Bird's eye view, focus skin
R = logical(G$m)
R[1755] = TRUE
descendants = get_descendants(R, G)
to_plot = descendants & ancestors

vertex_shapes = "square"
frame_colors = rep("black", G_anc$m)
frame_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "darkcyan"
frame_colors[!to_plot[ancestors]] = "lightgray"
fill_colors = rep("white", G_anc$m)
fill_colors[BH_rejections[ancestors]] = "gray45"
fill_colors[BH_outer_nodes[ancestors]] = "darkcyan"
fill_colors[BH_outer_nodes[ancestors] &! FBH_outer_nodes[ancestors]] = "white"
fill_colors[!to_plot[ancestors]] = "lightgray"
vertex_labels = rep("", G_anc$m)
#vertex_labels[BH_rejections[ancestors]] = invPerm(ord)
#vertex_labels = NULL
vertex_size = 5

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/output_ICD_tree_focus_skin.pdf"
pdf(file = plot_filename, width = 16, height = 10)
par(mar=c(0,0,0,0)) # make smaller margins
plot_tree(G_anc, vertex_shapes, frame_colors, fill_colors, vertex_labels, vertex_size)
dev.off()
command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)

###############################################################################

# Focus in on 1756: Diseases of Musculoskeletal System

R = logical(G$m)
R[1756] = TRUE
descendants = get_descendants(R, G)
to_plot = descendants & ancestors
G_to_plot = subset_graph(to_plot,G)

m = G_to_plot$m
children = G_to_plot$C
parents = G_to_plot$Pa
edges = c()
for(j in 1:m){
  for(child in children[[j]]){
    edges = c(edges, j, child)
  }
}
g<-graph(edges, n=m, directed=TRUE)
roots = which(sapply(G_to_plot$Pa, length) == 0)


vertex_shapes = "square"
frame_colors = rep("black", G_to_plot$m)
frame_colors[BH_outer_nodes[to_plot] &! FBH_outer_nodes[to_plot]] = "darkcyan"
fill_colors = rep("white", G_to_plot$m)
fill_colors[BH_rejections[to_plot]] = "gray45"
fill_colors[BH_outer_nodes[to_plot]] = "darkcyan"
fill_colors[BH_outer_nodes[to_plot] &! FBH_outer_nodes[to_plot]] = "white"
#vertex_labels = tree$meaning[to_plot]
vertex_labels = c("Infectious arthropathies", 
                  "Inflammatory polyarthropathies",
                  "Spondylopathies",
                  "Disorders of synovium and tendon",
                  "Diseases of the musculoskeletal system and connective tissue",
                  "Other rheumatoid arthritis",
                  "Rheumatoid arthritis, unspecified",
                  "RA (Multiple sites)",
                  "RA (Shoulder region)",
                  "RA (Hand)",
                  "Ankylosing spondylitis",
                  "AS (Site unspecified)",
                  "Synovitis and tenosynovitis",
                  "Other synovitis and tenosynovitis")
vertex_sizes = 50*strwidth(vertex_labels)

#vertex_labels = rep("", G_to_plot$m)
#vertex_labels[BH_rejections[ancestors]] = invPerm(ord)
#vertex_labels = NULL
vertex_size = 5

layout = layout_as_tree(g, root = roots)
#layout[c(1,2),] = layout[c(2,1),]

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/output_ICD_tree_AS.pdf"
pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
#     vertex.size = 15,
#     vertex.size2 = 5,
     vertex.size = vertex_sizes, 
     vertex.size2 = 15, 
     edge.width = 2, arrow.size = 2, 
     #     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
     #     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black",
     vertex.frame.color = frame_colors, layout = layout, asp = 0)
dev.off()
command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)

###############################################################################

# Focus in on 1755: Diseases of skin and subcutaneous tissue

R = logical(G$m)
R[1755] = TRUE
descendants = get_descendants(R, G)
to_plot = descendants & ancestors
G_to_plot = subset_graph(to_plot,G)

m = G_to_plot$m
children = G_to_plot$C
parents = G_to_plot$Pa
edges = c()
for(j in 1:m){
  for(child in children[[j]]){
    edges = c(edges, j, child)
  }
}
g<-graph(edges, n=m, directed=TRUE)
roots = which(sapply(G_to_plot$Pa, length) == 0)


vertex_shapes = "square"
frame_colors = rep("black", G_to_plot$m)
frame_colors[BH_outer_nodes[to_plot] &! FBH_outer_nodes[to_plot]] = "darkcyan"
fill_colors = rep("white", G_to_plot$m)
fill_colors[BH_rejections[to_plot]] = "gray45"
fill_colors[BH_outer_nodes[to_plot]] = "darkcyan"
fill_colors[BH_outer_nodes[to_plot] &! FBH_outer_nodes[to_plot]] = "white"
#vertex_labels = tree$meaning[to_plot]
vertex_labels = c("Infections of the skin and subcutaneous tissue",
                  "Papulosquamous disorders",
                  "Urticaria and erythema",
                  "Diseases of the skin and subcutaneous tissue",
                  "Cellulitis",
                  "Cellulitis of other sites",
                  "Psoriasis",
                  "Arthropathic psoriasis",
                  "Urticaria")
vertex_sizes = 60*strwidth(vertex_labels)

#vertex_labels = rep("", G_to_plot$m)
#vertex_labels[BH_rejections[ancestors]] = invPerm(ord)
#vertex_labels = NULL
vertex_size = 5

layout = layout_as_tree(g, root = roots)
#layout[c(1,2),] = layout[c(2,1),]

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/output_ICD_tree_skin.pdf"
pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
     #     vertex.size = 15,
     #     vertex.size2 = 5,
     vertex.size = vertex_sizes, 
     vertex.size2 = 15, 
     edge.width = 2, arrow.size = 2, 
     #     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
     #     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black",
     vertex.frame.color = frame_colors, layout = layout, asp = 0)
dev.off()
command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)



