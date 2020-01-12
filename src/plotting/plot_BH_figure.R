###########################################################
# CREATE PLOT SHOWING RESULTS OF LEAF BH AND BH EXPERIMENTS
###########################################################

# setup
machine = "local"
source("setup.R")

# names of PheWAS experiments
experiment_names = c("PheWAS_intermediate_leaf", "PheWAS_intermediate", "GO_Simes")

# read all results from file
results_all = vector("list", length(experiment_names))
for(experiment_idx in 1:length(experiment_names)){
  experiment_name = experiment_names[experiment_idx]
  results_dir = sprintf("%s/results/%s", base_dir, experiment_name)
  results_files = list.files(results_dir)
  results = vector("list", length(results_files))
  for(index in 1:length(results_files)){
    results[[index]] = read_tsv(sprintf("%s/%s", results_dir, results_files[[index]]), 
                                col_types = "ccddi", 
                                comment = "#")
  }
  df_results = do.call("rbind", results)
  if(experiment_name == "GO_Simes"){
    df_results$experiment = "GO"    
  } else{
    df_results$experiment = "ICD"
  }

  results_all[[experiment_idx]] = df_results
}

# massage all results into a data frame
df_results_all = do.call("rbind", results_all)

df_results_all = df_results_all %>% 
  filter(method %in% c("Focused_BH_original", "BH", "Leaf_BH_1", "Leaf_BH_2", "Leaf_BH_3")) %>%
  group_by(experiment, metric, method, signal_strength) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(n())) %>% 
  ungroup() %>%
  group_by(experiment, metric) %>%
  mutate(value_mean = ifelse(metric == "power", value_mean/max(value_mean), value_mean)) %>%
  ungroup() %>% 
  mutate(method = factor(method, 
                         levels = c("Focused_BH_original", 
                                    "BH",
                                    "Leaf_BH_1",
                                    "Leaf_BH_2",
                                    "Leaf_BH_3"),
                         labels = c("Focused BH",
                                    "BH (all nodes)",
                                    "BH (height 1)",
                                    "BH (height 2)",
                                    "BH (height 3)")),
         metric = factor(metric, 
                         levels = c("power", "fdp"), 
                         labels = c("Power", "False discovery rate")),
         experiment = factor(experiment, levels = c("ICD", "GO")))
        

# plot results and save to file
p = df_results_all %>% filter(signal_strength <= 6) %>% arrange(desc(method)) %>%
  ggplot(aes(x = signal_strength, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  geom_hline(data = tibble(metric = factor("False discovery rate", levels = c("Power", "False discovery rate")), q = 0.1),
             aes(yintercept = q), linetype = "dashed") +
  geom_errorbar(data = df_results_all %>% filter(metric == "False discovery rate", signal_strength <= 6),
                aes(ymin = value_mean - 2*value_se, ymax = value_mean + 2*value_se), 
                width = 0.2) + 
  facet_grid(metric ~ experiment, scales = "free_y") +
  scale_colour_manual(values = c("dodgerblue", "firebrick1", "gold2", "goldenrod3", "goldenrod4")) + 
  guides(colour = guide_legend(nrow = 2, byrow = FALSE)) + 
  # scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dashed")) + 
  theme_bw() + xlab("Signal strength") + ylab("Number of Power") +
  theme(legend.title = element_blank(), 
        legend.position = "bottom", axis.title.y = element_blank())
plot(p)
plot_filename = sprintf("%s/BH_experiments.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 5, height = 5)
