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
  

df %>% 
  group_by(depth, nonnull) %>%
  ggplot(aes(x = depth, group = nonnull, fill = nonnull)) + 
  geom_bar() + theme_bw() + scale_y_log10()

df_all %>% 
  mutate(experiment = factor(experiment, levels = experiment_names, labels = c("Clustered", "Intermediate", "Dispersed"))) %>%
  group_by(experiment, height, nonnull) %>%
  summarise(count = n()) %>% 
  ggplot(aes(x = height, y = log10(1+count), group = nonnull, fill = nonnull)) + 
  geom_col() + facet_wrap(experiment ~ .) + theme_bw() + 
  # scale_y_continuous(trans = "log1p", breaks = c(1,10,100,1000, 10000, 100000)) + 
  coord_flip()

df %>% 
  group_by(height, nonnull) %>%
  ggplot(aes(y = height, group = nonnull, fill = nonnull)) + 
  geom_bar() + theme_bw() + scale_x_log10()

df %>% 
  group_by(height) %>%
  ggplot(aes(x = height, y = items_per_node)) + 
  geom_violin() + theme_bw() + scale_y_log10() + coord_flip()

