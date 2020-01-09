#################################################################
# CREATE TABLES AND PLOTS SHOWING RESULTS OF PHEWAS DATA ANALYSIS
#################################################################

# setup
machine = "local"
source("setup.R")

# filenames
rejections_filename = sprintf("%s/data/processed/biobank/PheWAS_rejections.tsv", base_dir)
graph_filename = sprintf("%s/data/processed/biobank/G_at_least_50_cases.Rda", base_dir)
HES_affected_file = sprintf("%s/data/raw/biobank/HES_affected.Rda", base_dir)
tree_file = sprintf("%s/data/raw/biobank/coding19.tsv", base_dir)
pvals_filename = sprintf("%s/data/processed/biobank/pvals_HES_B_2705_m50.tsv", base_dir)

# load data
df_rejections = read_tsv(rejections_filename, col_types = "ccclc")
load(graph_filename)
G = G_reduced
load(HES_affected_file)
tree = read_tsv(tree_file, col_types = "cciic")
P = read_tsv(pvals_filename, col_types = "d") %>% pull()

############# VENN DIAGRAM ###################

# read rejections
FBH = df_rejections %>% 
  filter(rejected, 
         filter_status == "after_filter", 
         method == "Focused_BH") %>% 
  pull(node) 
SH = df_rejections %>% 
  filter(rejected, 
         filter_status == "after_filter", 
         method == "Structured_Holm") %>% 
  pull(node) 
LBH = df_rejections %>% 
  filter(rejected, 
         filter_status == "after_filter", 
         method == "Leaf_BH") %>% 
  pull(node) 

# create the plot
venn_plot_filename =  sprintf("%s/Venn_PheWAS.png", figures_dir)
venn.plot <- venn.diagram(
  x = list(
    LBH = LBH,
    FBH = FBH,
    SH = SH
  ),
  filename = venn_plot_filename,
  col = "transparent",
  fill = c("gold2", "dodgerblue", "darkorange2"),
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  # euler.d = TRUE,
  scaled = TRUE,
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = c(0.07, 0.07, 0.02),
  cat.fontfamily = "serif",
  margin = 0.3,
  imagetype = "png"
);

# trim whitespace
system(sprintf("convert %s -trim %s", venn_plot_filename, venn_plot_filename))

############# TABLE OF DISCOVERIES BY DISEASE CLASS ###################

# add various structural information to table
P = P[P < 1]
names(P) = G$node_names
node_depths = get_depths(G$Pa)
names(node_depths) = G$node_names
node_heights = get_heights(G$C)
names(node_heights) = G$node_names
nodes = unique(df_rejections$node)
num_affected = sapply(nodes, function(node)(length(affected[[as.character(tree$node_id[tree$coding == node])]])))
root_nodes = get_root_nodes(G$Pa)
names(root_nodes) = G$node_names
roots = which(node_depths == 1)
root_meanings = tree %>% filter(coding %in% names(roots)) %>% pull(meaning)
names(root_meanings) = tree %>% filter(coding %in% names(roots)) %>% pull(coding)
root_short_meanings = c("Infectious", "Neoplasms", "Blood", "Endocrine", 
                        "Circulatory", "Psychological", "Nervous", "Eye",
                        "Ear", "Respiratory", "Digestive", "Skin", "Musculoskeletal",
                        "Genitourinary", "Injury", "Pregnancy", "Congenital",
                        "Clinical symptoms", "Health services")
names(root_short_meanings) = tree %>% filter(coding %in% names(roots)) %>% pull(coding)

df_rejections = df_rejections %>% 
  mutate(depth = node_depths[node], 
         height = node_heights[node],
         affected = num_affected[node],
         category = root_short_meanings[G$node_names[root_nodes[node]]],
         P = P[node])

df_rejections %>% 
  filter(filter_status == "after_filter", category %in% categories_with_rejections) %>% 
  group_by(method, category) %>% 
  summarise(rejections_per_category = sum(rejected)) %>%
  ungroup() %>%
  spread(method, rejections_per_category) %>%
  adorn_totals("row") %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
        col.names = c("",
                      "BH",
                      "FBH",
                      "LBH", 
                      "SH",
                      "Y")) %>% 
  row_spec(11, hline_after = TRUE) %>%
  save_kable(sprintf("%s/PheWAS_table.pdf", figures_dir))

############# TREE PLOT OF MUSCULOSKELETAL DISEASES ###################

FBH_rejections = df_rejections %>% filter(filter_status == "after_filter", 
                                          method == "Focused_BH", 
                                          category == "Musculoskeletal") %>% pull(rejected)
SH_rejections = df_rejections %>% filter(filter_status == "after_filter", 
                                          method == "Structured_Holm", 
                                          category == "Musculoskeletal") %>% pull(rejected)
Leaf_BH_rejections = df_rejections %>% filter(filter_status == "after_filter", 
                                         method == "Leaf_BH", 
                                         category == "Musculoskeletal") %>% pull(rejected)
musculoskeletal_nodes = root_short_meanings[G$node_names[root_nodes]] == "Musculoskeletal"
G_musculo = subset_graph(!musculoskeletal_nodes, G)
leaf_nodes = sapply(G_musculo$C, length) == 0
FBH_relatives = get_ancestors(FBH_rejections, G_musculo) | get_descendants(FBH_rejections, G_musculo$C)
leaf_nodes_to_plot = logical(G_musculo$m)
leaf_nodes_to_plot[FBH_relatives & leaf_nodes] = TRUE
other_leaf_nodes = !FBH_relatives & leaf_nodes & get_descendants(FBH_relatives & node_depths[musculoskeletal_nodes] > 1,G_musculo$C)
downsampling_factor = 3
leaf_nodes_to_plot[withSeed(sample(which(other_leaf_nodes), floor(sum(other_leaf_nodes)/downsampling_factor)),1)] = TRUE
nodes_to_plot = get_ancestors(leaf_nodes_to_plot, G_musculo)

G_to_plot = subset_graph(!nodes_to_plot, G_musculo)

edges = c()
for(j in 1:G_to_plot$m){
  for(child in G_to_plot$C[[j]]){
    edges = c(edges, j, child)
  }
}
g<-graph(edges, n=G_to_plot$m, directed=TRUE)

tree %>% filter(coding %in% G_to_plot$node_names) %>% pull(meaning)

node_labels = rep("", G_to_plot$m)
node_labels[FBH_rejections[nodes_to_plot]] = rank(P[musculoskeletal_nodes][FBH_rejections])
V(g)$node_label <- node_labels
V(g)$`Rejected by` = rep("none", G_to_plot$m)
V(g)$`Rejected by`[FBH_rejections[nodes_to_plot] & SH_rejections[nodes_to_plot] & Leaf_BH_rejections[nodes_to_plot]] = "FBH, SH, LBH"
V(g)$`Rejected by`[FBH_rejections[nodes_to_plot] & SH_rejections[nodes_to_plot] & !Leaf_BH_rejections[nodes_to_plot]] = "FBH, SH"
V(g)$`Rejected by`[FBH_rejections[nodes_to_plot] & !SH_rejections[nodes_to_plot] & Leaf_BH_rejections[nodes_to_plot]] = "FBH, LBH"
V(g)$`Rejected by`[FBH_rejections[nodes_to_plot] & !SH_rejections[nodes_to_plot] & !Leaf_BH_rejections[nodes_to_plot]] = "FBH"
V(g)$`Disease cases` = num_affected[musculoskeletal_nodes][nodes_to_plot]
V(g)$`p-value` = 1/pmin(1e-2, pmax(1e-6, P[musculoskeletal_nodes][nodes_to_plot]))

p = ggraph(g,layout='tree')+
  geom_edge_link(arrow=arrow(length=unit(2,'mm')),end_cap=circle(3,'mm'))+
  geom_node_point(aes(size = `p-value`, colour = `Rejected by`)) + 
  geom_node_text(aes(label=node_label), nudge_x = 0, nudge_y = 0) +
  scale_radius(trans = "log10", 
               range = c(4,8),
               breaks = c(1e2, 1e3, 1e4, 1e5, 1e6),
               labels = c("<1", "<1e-2", "<1e-3", "<1e-4", "<1e-5")) + 
  scale_colour_manual(values = c("seagreen3","gold2","darkorange2", "dodgerblue", "grey50"), 
                      breaks = c("FBH, SH, LBH", "FBH, SH", "FBH, LBH", "FBH")) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  theme_void() + ggtitle("Musculoskeletal diseases\nassociated with HLA allele B*27:05") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "left")
plot(p)
plot_filename = sprintf("%s/musculoskeletal_graph.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 5, height = 4)

############# TABLE OF DISCOVERED MUSCULOSKELETAL DISEASES ###################

annotation_table = df_rejections %>% filter(filter_status == "after_filter", 
                                            method == "Focused_BH", 
                                            category == "Musculoskeletal",
                                            rejected) %>% 
  select(meaning, affected, P)
annotation_table$meaning = c(linebreak("Infectious arthropathies"), 
                             linebreak("Rheumatoid arthritis\n(multiple sites)"), 
                             linebreak("Rheumatoid arthritis\n(shoulder)"), 
                             linebreak("Rheumatoid arthritis\n(hand)"), 
                             "Ankylosing spondylitis", 
                             linebreak("Ankylosing spondylitis\n(unspecified)"), 
                             "Synovitis/tenosynovitis")
annotation_table = annotation_table %>% arrange(P) %>%
  mutate(P = format(P, scientific = TRUE, digits = 3),
         index = row_number()) %>%
  rename(Disease = meaning, Cases = affected, `p-value` = P) %>% 
  select(index, Disease, Cases, `p-value`) 
  
annotation_table %>% kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
                           col.names = c("",
                                         "Disease",
                                         "Cases",
                                         "p-value")) %>% 
  save_kable(sprintf("%s/musculoskeletal_table.pdf", figures_dir))
