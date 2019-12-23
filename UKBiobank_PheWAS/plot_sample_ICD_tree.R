rm(list = ls())
data_dir = "/home/ekatsevi/Documents/large_file_storage/filtered_BH/biobank"
setwd("/home/ekatsevi/Dropbox/Stanford/Research/Projects/HierTest/code/UKBiobank_PheWAS")
source("aux_tree.R")
library(readr)
library(Matrix)
options(tibble.print_max = Inf) 

# read tree
tree_file = sprintf("%s/coding19.tsv", data_dir)
tree = read_tsv(tree_file)
m = nrow(tree)

# read graph structure
graph_file = sprintf("%s/graph.Rda", data_dir)
load(graph_file)

to_plot = logical(m)
to_plot[1104] = TRUE # Symptoms and signs involving the circulatory and respiratory systems (level 1)
to_plot[12291] = TRUE # Abnormalities of breathing (level 2)
to_plot[12271] = TRUE # Abnormalities of heart beat (level 2)
to_plot[12301] = TRUE # Pain in throat and chest (level 2)
to_plot[12272] = TRUE # Tachycardia (level 3)
to_plot[12274] = TRUE # Palpitations (level 3)
to_plot[12294] = TRUE # Wheezing (level 3)
to_plot[12297] = TRUE # Mouth breathing (level 3)
# to_plot[12298] = TRUE # Hiccough (level 3)
to_plot[12302] = TRUE # Pain in throat (level 3)
to_plot[12302] = TRUE # Pain in throat (level 3)
to_plot[12306] = TRUE # Chest pain, unspecified (level 3)

G_to_plot = subset_graph(to_plot,G)

vertex_shapes = "rectangle"
frame_colors = NULL
#fill_colors[6] = "red"
vertex_labels = tree$meaning[to_plot]
#vertex_labels = NULL
vertex_size = 5

library(igraph)
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


vertex_labels[3] = "R00.0 Tachycardia"
vertex_labels[10] = "R07.4 Chest pain"
vertex_sizes = 45*strwidth(vertex_labels)
fill_colors = rep("white", G_to_plot$m)
fill_colors[7] = "darkcyan"

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/sample_ICD_tree.pdf"
pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
     vertex.size = vertex_sizes, vertex.size2 = 15, edge.width = 2, arrow.size = 2, vertex.frame.width = 5,
#     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
#     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black",
     vertex.frame.color = frame_colors, layout = layout_as_tree(g, root = roots), asp = 0)
dev.off()

command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)

par(mar=c(5.1, 4.1, 4.1, 2.1)) # reset margins

#########################################

fill_colors = rep("White", G_to_plot$m)
fill_colors[c(8)] = "firebrick1"
fill_colors[c(1,2,3)] = "darkcyan"

vertex_frame_colors = rep("firebrick1", G_to_plot$m)
vertex_frame_colors[c(1,2,3, 4)] = "darkcyan"

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/outer_nodes_problem_1.pdf"

pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
     vertex.size = vertex_sizes, vertex.size2 = 15, edge.width = 2, arrow.size = 2, vertex.frame.width = 5,
     #     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
     #     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black", vertex.frame.color = vertex_frame_colors,
     vertex.frame.color = frame_colors, layout = layout_as_tree(g, root = roots), asp = 0)
dev.off()

command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)

#########################################################

fill_colors = rep("White", G_to_plot$m)
fill_colors[8] = "firebrick1"
fill_colors[3] = "darkcyan"
fill_colors[c(1,2)] = "Gray"

vertex_frame_colors = rep("firebrick1", G_to_plot$m)
vertex_frame_colors[c(1,2,3, 4)] = "darkcyan"

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/outer_nodes_problem_2.pdf"

pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
     vertex.size = vertex_sizes, vertex.size2 = 15, edge.width = 2, arrow.size = 2, vertex.frame.width = 5,
     #     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
     #     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black", vertex.frame.color = vertex_frame_colors,
     vertex.frame.color = frame_colors, layout = layout_as_tree(g, root = roots), asp = 0)
dev.off()

command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)


#########################################################

fill_colors = rep("White", G_to_plot$m)
fill_colors[c(1,2,3,8)] = "darkcyan"
vertex_frame_colors = rep("black", G_to_plot$m)

#vertex_frame_colors = rep("darkcyan", G_to_plot$m)
#vertex_frame_colors[c(1,2,3, 4)] = "darkcyan"

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/outer_nodes_illustration_1.pdf"

pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
     vertex.size = vertex_sizes, vertex.size2 = 15, edge.width = 2, arrow.size = 2, vertex.frame.width = 5,
     #     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
     #     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black", vertex.frame.color = vertex_frame_colors,
     vertex.frame.color = frame_colors, layout = layout_as_tree(g, root = roots), asp = 0)
dev.off()

command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)


#########################################################

fill_colors = rep("White", G_to_plot$m)
fill_colors[c(3,8)] = "darkcyan"
fill_colors[c(1,2)] = "gray"
vertex_frame_colors = rep("black", G_to_plot$m)

#vertex_frame_colors = rep("darkcyan", G_to_plot$m)
#vertex_frame_colors[c(1,2,3, 4)] = "darkcyan"

plot_filename = "/home/ekatsevi/Dropbox/Job_Search/Interview_Preparation/Job_Talk/figures/outer_nodes_illustration_2.pdf"

pdf(file = plot_filename, width = 16, height = 10)
#pdf(file = plot_filename, width = 12, height = 8)
par(mar=c(0,0,0,0)) # make smaller margins
plot(g, vertex.shape = "rectangle", vertex.label = vertex_labels,
     vertex.size = vertex_sizes, vertex.size2 = 15, edge.width = 2, arrow.size = 2, vertex.frame.width = 5,
     #     vertex.size=(strwidth(tree$meaning[to_plot]) + strwidth("oo")) * 10,
     #     vertex.size2=strheight("I") * 2 * 10)
     vertex.color = fill_colors, vertex.label.color = "black", vertex.frame.color = vertex_frame_colors,
     vertex.frame.color = frame_colors, layout = layout_as_tree(g, root = roots), asp = 0)
dev.off()

command = sprintf("pdfcrop %s %s", plot_filename, plot_filename)
system(command)





