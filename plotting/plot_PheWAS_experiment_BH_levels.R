###################################################
# CREATE PLOT SHOWING RESULTS OF LEAF BH EXPERIMENT
###################################################

# setup
machine = "local"
source("setup.R")

# name of experiment
experiment_name = "PheWAS_intermediate_leaf"

# read results from file
results_dir = sprintf("%s/results/%s", base_dir, experiment_name)
results_files = list.files(results_dir)
results = vector("list", length(results_files))
for(index in 1:length(results_files)){
  results[[index]] = read_tsv(sprintf("%s/%s", results_dir, results_files[[index]]), 
                              col_types = "ccddi", 
                              comment = "#")
}

# massage results into data frame, rename some entries
df_results = do.call("rbind", results)
df_results = df_results %>% mutate(method = factor(method, 
                                                   levels = c(
                                                              "Focused_BH_original", 
                                                              "Leaf_BH_1",
                                                              "Leaf_BH_2",
                                                              "Leaf_BH_3"),
                                                   labels = c(
                                                              "Focused BH",
                                                              "BH at height 1",
                                                              # "Focused BH (permutation)",
                                                              "BH at height 2",
                                                              "BH at height 3")),
                                   metric = factor(metric, 
                                                   levels = c("fdp", "power"), 
                                                   labels = c("False discovery rate", "Number of filtered discoveries"))) %>%
  group_by(method, metric, signal_strength) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(n()))

# plot results and save to file
p = df_results %>% filter(metric == "Number of filtered discoveries", signal_strength <= 6) %>% arrange(desc(method)) %>%
  ggplot(aes(x = signal_strength, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  scale_colour_manual(values = c("dodgerblue", "gold2", "goldenrod3", "goldenrod4")) + 
  # scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dashed")) + 
  theme_bw() + xlab("Signal strength") + ylab("Number of true discoveries") +
  theme(legend.title = element_blank(), 
        legend.position = c(0.225,0.675))
plot(p)
plot_filename = sprintf("%s/Leaf_BH.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 3.5, height = 2.5)