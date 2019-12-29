figures_dir = "/home/ekatsevi/Dropbox/Research/Projects/HierTest/manuscript/Reresubmission/figures"
# experiment_name = "PheWAS_clustered"
experiment_name = "PheWAS_dispersed"
results_dir = sprintf("%s/results/%s", base_dir, experiment_name)
results_files = list.files(results_dir)
results = vector("list", length(results_files))
for(index in 1:length(results_files)){
  results[[index]] = read_tsv(sprintf("%s/%s", results_dir, results_files[[index]]), 
                              col_types = "ccddi", 
                              comment = "#")
}
df_results = do.call("rbind", results)
df_results = df_results %>% mutate(method = factor(method, 
                                      levels = c("BH", 
                                                 "Focused_BH_original", 
                                                 "Leaf_BH",
                                                 # "Focused_BH_permutation", 
                                                 "Structured_Holm",
                                                 "Yekutieli"),
                                      labels = c("BH", 
                                                 "Focused BH", 
                                                 "Leaf BH",
                                                 # "Focused BH (permutation)", 
                                                 "Structured Holm",
                                                 "Yekutieli")),
                      metric = factor(metric, 
                                      levels = c("fdp", "power"), 
                                      labels = c("False discovery rate", "Number of filtered discoveries"))) %>%
  group_by(method, metric, signal_strength) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(n()))
p = df_results %>%
  ggplot(aes(x = signal_strength, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  scale_colour_manual(values = c("firebrick1", "dodgerblue", "blue", "darkgoldenrod", "purple")) + 
  geom_hline(data = tibble(metric = "False discovery rate", q = 0.1), 
             aes(yintercept = q), linetype = "dashed") + 
  geom_errorbar(data = df_results %>% filter(metric == "False discovery rate"),
                aes(ymin = value_mean - value_se, ymax = value_mean + value_se), 
                width = 0.2) + 
  facet_wrap(metric ~ ., scales = "free") + theme_bw() + xlab("Signal strength") +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = "bottom")
plot(p)
# plot_filename = sprintf("%s/REVIGO_experiment.pdf", figures_dir)
# ggsave(filename = plot_filename, plot = p, device = "pdf", width = 6, height = 3.5)