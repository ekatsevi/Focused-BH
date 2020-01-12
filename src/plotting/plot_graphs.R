rm(list = ls())
machine = "local"
source("setup.R")

experiment_names = c("PheWAS_clustered", "PheWAS_intermediate", "PheWAS_dispersed")
list_all = vector("list", length(experiment_names))
names(list_all) = experiment_names
for(experiment_name in experiment_names){
  input_filename = sprintf("simulations/input_files/input_file_%s.R", experiment_name)
  input_mode = "experiment"
  experiment_index = 1
  source(input_filename, local = TRUE)
  
  depths = get_depths(G$Pa)
  heights = get_heights(G$C)
  
  df = tibble(depth = as.factor(depths), height = as.factor(heights), items_per_node = items_per_node, nonnull = nonnull_nodes)
  df$experiment = experiment_name
  list_all[[experiment_name]] = df
}
df_all = do.call("rbind", list_all)
  
p_nodes = df_all %>% 
  mutate(experiment = factor(experiment, levels = experiment_names, labels = c("Clustered", "Intermediate", "Dispersed"))) %>%
  group_by(experiment, height) %>%
  summarise(nonnull_count = sum(nonnull), total_count = n()) %>% 
  ungroup() %>%
  spread(experiment, nonnull_count) %>%
  gather(counting_what, count, -height) %>%
  mutate(counting_what = factor(counting_what, 
                                levels = c("Clustered", "Intermediate", "Dispersed", "total_count"),
                                labels = c("Non-nulls (clustered)", "Non-nulls (intermediate)", "Non-nulls (dispersed)", "Total nodes"))) %>%
  ggplot(aes(x = height, y = count, group = counting_what, fill = counting_what)) + 
  geom_col(position = position_dodge2(reverse = TRUE)) + 
  scale_fill_manual(values = c("cyan3", "deepskyblue4",  "navy", "dimgray")) + 
  theme_bw() + xlab("Height in graph") + ylab("Number of nodes") + 
  scale_y_continuous(trans = "log1p", breaks = c(1,10,100,1000,10000)) + 
  coord_flip() + theme(legend.title = element_blank(), legend.position = c(0.99,0.99), legend.justification = c(1,1),#legend.position = c(0.725, 0.825),
                       panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
plot(p_nodes)

p_items = df_all %>% 
  filter(experiment == "PheWAS_clustered") %>%
  ggplot(aes(x = height, y = items_per_node)) + 
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_log10() + coord_flip()  + 
  theme_bw() + ylab("Number of items per node") + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title.y = element_blank(), 
        # axis.text.y = element_blank()
        )
plot(p_items)

p_combined = grid.arrange(p_nodes, p_items, ncol = 2)

plot_filename = sprintf("%s/ICD_simulation_graph_info.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p_combined, device = "pdf", width = 8, height = 3)

###################################


experiment_name = "GO_Simes"
input_filename = sprintf("simulations/input_files/input_file_%s.R", experiment_name)
input_mode = "experiment"
experiment_index = 1
source(input_filename, local = TRUE)

depths = get_depths(G$Pa)
heights = get_heights(G$C)

df = tibble(depth = as.factor(depths), height = as.factor(heights), items_per_node = items_per_node, nonnull = nonnull_nodes)

p_nodes = df %>% 
  group_by(height) %>%
  summarise(nonnull_count = sum(nonnull), total_count = n()) %>% 
  ungroup() %>%
  gather(counting_what, count, -height) %>%
  mutate(counting_what = factor(counting_what, 
                                levels = c("nonnull_count", "total_count"),
                                labels = c("Non-null nodes", "Total nodes"))) %>%
  ggplot(aes(x = height, y = count, group = counting_what, fill = counting_what)) + 
  geom_col(position = position_dodge2(reverse = TRUE)) + 
  scale_fill_manual(values = c("navy", "dimgray")) + 
  theme_bw() + xlab("Height in graph") + ylab("Number of nodes") + 
  scale_y_continuous(trans = "log1p", breaks = c(1,10,100,1000)) + 
  coord_flip() + theme(legend.title = element_blank(), legend.position = c(0.99,0.99), legend.justification = c(1,1),#legend.position = c(0.725, 0.825),
                       panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
plot(p_nodes)

p_items = df %>% 
  ggplot(aes(x = height, y = items_per_node)) + 
  # geom_violin(scale = "width") + 
  geom_boxplot(outlier.shape = NA) + 
  scale_y_log10() + 
  coord_flip()  + 
  theme_bw() + ylab("Number of items per node") + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title.y = element_blank(), 
        # axis.text.y = element_blank()
        )
plot(p_items)

p_combined = grid.arrange(p_nodes, p_items, ncol = 2)

plot_filename = sprintf("%s/REVIGO_experiment_graph_info.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p_combined, device = "pdf", width = 8, height = 3)

###################################

# preprocessed GO graph
reduced_graph_file = sprintf("%s/data/processed/GO/GO_reduced_graph.Rda", base_dir)
load(reduced_graph_file)

# create adjacency matrix between genes and GO terms
num_items = G$num_items
m = G$m
adj_matrix = matrix(FALSE, num_items, m)
for(node in 1:m){
  genes_in_node = G$sets[[node]]
  adj_matrix[genes_in_node, node] = TRUE
}
items_per_node = sapply(G$sets, length)

depths = get_depths(G$Pa)
heights = get_heights(G$C)

df = tibble(depth = factor(depths, levels = 1:max(depths)), height = as.factor(heights), items_per_node = items_per_node)


p_nodes = df %>%
  group_by(height) %>%
  summarise(total_count = n()) %>%
  ggplot(aes(x = height, y = total_count)) +
  geom_col(fill = "dimgray") +
  theme_bw() + xlab("Height in graph") + ylab("Number of nodes") +
  scale_y_continuous(trans = "log1p", breaks = c(1,10,100,1000)) +
  coord_flip() + theme(legend.title = element_blank(), legend.position = c(0.99,0.99), legend.justification = c(1,1),#legend.position = c(0.725, 0.825),
                       panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
plot(p_nodes)

p_items = df %>% 
  ggplot(aes(x = height, y = items_per_node)) + 
  # geom_violin(scale = "width") +
  # geom_point(data = df %>% filter(height %in% c(17, 18))) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_y_log10() + 
  coord_flip()  + 
  # geom_hline(yintercept = 100, linetype = "dashed", colour = "red") + 
  theme_bw() + ylab("Number of items per node") + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),        
        # axis.text.y = element_blank(),
        axis.title.y = element_blank())
plot(p_items)

p_combined = grid.arrange(p_nodes, p_items, ncol = 2)

plot_filename = sprintf("%s/REVIGO_analysis_graph_info.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p_combined, device = "pdf", width = 8, height = 4)


###################################

# preprocessed ICD graph
min_affected = 50
reduced_graph_filename = sprintf("%s/data/processed/biobank/G_at_least_%d_cases.Rda", base_dir, min_affected)
load(reduced_graph_filename)
G = G_reduced

# create adjacency matrix between genes and GO terms
num_items = G$num_items
m = G$m
adj_matrix = matrix(FALSE, num_items, m)
for(node in 1:m){
  genes_in_node = G$sets[[node]]
  adj_matrix[genes_in_node, node] = TRUE
}
items_per_node = sapply(G$sets, length)

depths = get_depths(G$Pa)
heights = get_heights(G$C)

df = tibble(depth = as.factor(depths), height = as.factor(heights), items_per_node = items_per_node)

p_nodes = df %>%
  group_by(height) %>%
  summarise(total_count = n()) %>%
  ggplot(aes(x = height, y = total_count)) +
  geom_col(fill = "dimgray") +
  theme_bw() + xlab("Height in graph") + ylab("Number of nodes") +
  scale_y_continuous(trans = "log1p", breaks = c(1,10,100,1000)) +
  coord_flip() + theme(legend.title = element_blank(), legend.position = c(0.99,0.99), legend.justification = c(1,1),#legend.position = c(0.725, 0.825),
                       panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
plot(p_nodes)

p_items = df %>% 
  ggplot(aes(x = height, y = items_per_node)) + 
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_log10() + 
  coord_flip()  + 
  theme_bw() + ylab("Number of items per node") + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title.y = element_blank(), 
        # axis.text.y = element_blank()
        )
plot(p_items)

p_combined = grid.arrange(p_nodes, p_items, ncol = 2)

plot_filename = sprintf("%s/ICD_analysis_graph_info.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p_combined, device = "pdf", width = 8, height = 2)

