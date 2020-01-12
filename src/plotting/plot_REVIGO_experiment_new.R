###################################################
# CREATE PLOT SHOWING RESULTS OF REVIGO EXPERIMENT
###################################################

# setup 
machine = "local"
source("setup.R")

# name of experiment
experiment_name = "GO_Simes"

# read results from file
results_dir = sprintf("%s/results/%s", base_dir, experiment_name)
results_files = list.files(results_dir)
results = vector("list", length(results_files))
for(index in 1:length(results_files)){
  results[[index]] = read_tsv(sprintf("%s/%s", results_dir, results_files[[index]]), 
                              col_types = "ccddi", 
                              comment = "#")
}

# massage results into a data frame, rename some column entries
df_results = do.call("rbind", results)
df_results = df_results %>% 
  filter(method != "BH") %>%
  group_by(method, metric, signal_strength) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(n())) %>%
  ungroup() %>%
  group_by(metric) %>%
  mutate(value_mean = ifelse(metric == "power", value_mean/max(value_mean), value_mean)) %>%
  ungroup() %>%
  mutate(method = factor(method, 
                         levels = c("Focused_BH_original", 
                                    "Focused_BH_permutation",
                                    "Structured_Holm"),
                         labels = c("Focused BH", 
                                    "Focused BH (permutation)",
                                    "Structured Holm")),
         metric = factor(metric, 
                         levels = c("fdp", "power"),
                         labels = c("False discovery rate", "Power")))
  

# plot the results and save to file
p = df_results %>% filter(signal_strength <= 6) %>%
  ggplot(aes(x = signal_strength, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  scale_colour_manual(values = c("dodgerblue", "blue", "darkorange2")) + 
  geom_hline(data = tibble(metric = "False discovery rate", q = 0.1), 
             aes(yintercept = q), linetype = "dashed") + 
  geom_errorbar(data = df_results %>% filter(metric == "False discovery rate", signal_strength <= 6),
                aes(ymin = value_mean - 2*value_se, ymax = value_mean + 2*value_se), 
                width = 0.2) + 
  facet_wrap(metric ~ ., scales = "free") + theme_bw() + xlab("Signal strength") +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = "bottom")
plot(p)
plot_filename = sprintf("%s/REVIGO_experiment_new.pdf", figures_dir)
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 6, height = 3.5)