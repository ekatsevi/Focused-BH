#############################################################
# CREATE TABLES AND PLOT ILLUSTRATING GO/REVIGO DATA ANALYSIS
#############################################################

# setup
machine = "local"
source("setup.R")

# filenames
rejections_filename = sprintf("%s/data/processed/GO/REVIGO_rejections.tsv", base_dir)
graph_filename = sprintf("%s/data/processed/GO/GO_reduced_graph.Rda", base_dir)
GO_data_filename = sprintf("%s/data/raw/GO/GO.xls", base_dir)

# read data from file
df_rejections = read_tsv(rejections_filename, col_types = "cccl")
load(graph_filename)
m = G$m
ids = G$node_names
GO_data = read_excel(GO_data_filename)
P = rep(1, m)
significant_ids = match(GO_data$`GO Term`, ids)
P[significant_ids[!is.na(significant_ids)]] = GO_data$`P-value`[!is.na(significant_ids)]
names(P) = ids

############# TABLE OF NUMBERS OF DISCOVERIES BEFORE AND AFTER FILTERING ###################
df_rejections %>% 
  filter(rejected) %>% 
  group_by(filter_status, method) %>% 
  summarise(num_rejections = sum(rejected)) %>% 
  ungroup() %>%
  spread(method, num_rejections) %>%
  mutate(filter_status = factor(filter_status, 
                                levels = c("before_filter", "after_filter"), 
                                labels = c("Before filtering", "After filtering"))) %>%
  arrange(filter_status) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
                                        col.names = c("",
                                                      "BH",
                                                      "Focused BH",
                                                      "Structured Holm")) %>%
  save_kable(sprintf("%s/REVIGO_summary_table.pdf", figures_dir))
  
############# TREE PLOT ###################
term_descriptions = GO_data$Description
names(term_descriptions) = GO_data$`GO Term`
df_rejections = df_rejections %>% 
  mutate(P = format(P[node], scientific = TRUE, digits = 3), description = term_descriptions[node])

FBH_rejections = df_rejections %>% filter(filter_status == "after_filter", 
                                          method == "Focused_BH") %>% pull(rejected)
SH_rejections = df_rejections %>% filter(filter_status == "after_filter", 
                                         method == "Structured_Holm") %>% pull(rejected)
all_rejections = FBH_rejections | SH_rejections
leaf_nodes = sapply(G$C, length) == 0
rejections_relatives = get_ancestors(all_rejections, G) # | get_descendants(all_rejections, G$C)
nodes_to_plot = rejections_relatives

G_to_plot = subset_graph(!nodes_to_plot, G)

edges = c()
for(j in 1:G_to_plot$m){
  for(child in G_to_plot$C[[j]]){
    edges = c(edges, j, child)
  }
}
g<-graph(edges, n=G_to_plot$m, directed=TRUE)

node_labels = rep("", G_to_plot$m)
node_labels[all_rejections[nodes_to_plot]] = rank(P[all_rejections])
V(g)$node_label <- node_labels
V(g)$`Rejected by` = rep("none", G_to_plot$m)
V(g)$`Rejected by`[FBH_rejections[nodes_to_plot] & SH_rejections[nodes_to_plot]] = "FBH, SH"
V(g)$`Rejected by`[FBH_rejections[nodes_to_plot] & !SH_rejections[nodes_to_plot]] = "FBH"
V(g)$`Rejected by`[!FBH_rejections[nodes_to_plot] & SH_rejections[nodes_to_plot]] = "SH"
V(g)$`p-value` = 1/pmin(1e-2, pmax(1e-6, P[rejections_relatives]))

p = ggraph(g,layout='sugiyama')+
  geom_edge_link(arrow=arrow(length=unit(2,'mm')),end_cap=circle(3,'mm'))+
  geom_node_point(aes(size = `p-value`, colour = `Rejected by`)) + 
  geom_node_text(aes(label=node_label), nudge_x = 0, nudge_y = 0) +
  scale_colour_manual(values = c("seagreen3","dodgerblue", "grey50", "darkorange2"), 
                      breaks = c("FBH, SH", "FBH", "SH")) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  scale_radius(trans = "log10", 
                       range = c(4,8),
                       breaks = c(1e2, 1e3, 1e4, 1e5, 1e6),
                       labels = c("<1", "<1e-2", "<1e-3", "<1e-4", "<1e-5")) +
  theme_void() + ggtitle("Biological processes associated with \nbreast cancer outcome") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "left")
plot(p)

plot_filename = sprintf("%s/GO_graph.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 5, height = 4)

############# TABLE OF REJECTED GO TERMS ###################
annotation_table = tibble(`Biological process` = term_descriptions[match(ids[all_rejections], GO_data$`GO Term`)],
                          `p-value` = P[all_rejections]) %>% arrange(`p-value`)

annotation_table$`Biological process`[4] = linebreak("negative regulation of\nchromosome segregation")
annotation_table$`Biological process`[9] = linebreak("anaphase-promoting complex-dependent\ncatabolic process")
annotation_table$`Biological process`[12] = linebreak("cellular component\norganization or biogenesis")

annotation_table = annotation_table %>% mutate(index = row_number(), `p-value` = format(`p-value`, scientific = TRUE, digits = 3)) %>%
  select(index, `Biological process`, `p-value`)

annotation_table %>% kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
                           col.names = c("",
                                         "Biological process",
                                         "p-value")) %>% 
  save_kable(sprintf("%s/REVIGO_table.pdf", figures_dir))
