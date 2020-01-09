######################################################
# CREATE PLOT SHOWING RESULTS OF PHEWAS EXPERIMENT
######################################################

# setup
machine = "local"
source("setup.R")

# names of PheWAS experiments in different signal regimes
experiment_names = c("PheWAS_clustered", "PheWAS_intermediate", "PheWAS_dispersed")

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
  df_results$experiment = experiment_name
  results_all[[experiment_idx]] = df_results
}

# massage all results into a data frame
df_results_all = do.call("rbind", results_all)
df_results_all = df_results_all %>% filter(method != "Leaf_BH")
df_results_all = df_results_all %>% mutate(method = factor(method, 
                                                   levels = c("BH", 
                                                              "Focused_BH_original", 
                                                              "Structured_Holm",
                                                              "Yekutieli"),
                                                   labels = c("BH", 
                                                              "Focused BH", 
                                                              "Structured Holm",
                                                              "Yekutieli")),
                                   metric = factor(metric, 
                                                   levels = c("power", "fdp"), 
                                                   labels = c("True discoveries", "False discovery rate")),
                                   experiment = factor(experiment,
                                                       levels = experiment_names,
                                                       labels = c("Clustered", "Intermediate", "Dispersed"))) %>%
  group_by(method, metric, experiment, signal_strength) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(n()))

# plot results and save to file
p = df_results_all %>% filter(signal_strength <= 6) %>%
  ggplot(aes(x = signal_strength, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  scale_colour_manual(values = c("firebrick1", "dodgerblue", 
                                  "darkorange2",
                                 # "darkgoldenrod", 
                                 "purple")) + 
  geom_hline(data = tibble(metric = factor("False discovery rate", levels = c("True discoveries", "False discovery rate")), q = 0.1),
             aes(yintercept = q), linetype = "dashed") +
  geom_errorbar(data = df_results_all %>% filter(metric == "False discovery rate", signal_strength <= 6),
                aes(ymin = value_mean - 2*value_se, ymax = value_mean + 2*value_se), 
                width = 0.2) + 
  facet_grid(metric ~ experiment, scales = "free_y") + theme_bw() + xlab("Signal strength") +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = "bottom")
plot(p)
plot_filename = sprintf("%s/PheWAS_experiment.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 6, height = 4.25)