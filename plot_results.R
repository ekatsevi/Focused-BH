experiment_name = "GO_Simes"
results_dir = sprintf("%s/results/%s", base_dir, experiment_name)
results_files = list.files(results_dir)
results = vector("list", length(results_files))
for(index in 1:length(results_files)){
  results[[index]] = read_tsv(sprintf("%s/%s", results_dir, results_files[[index]]), col_types = "ccidd")
}
df_results = do.call("rbind", results)
df_results %>% group_by(method, metric, signal_strength) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = signal_strength, y = value_mean, group = method)) + 
  geom_line(aes(colour = method)) + 
  geom_errorbar(aes(ymin = value_mean - value_se, ymax = value_mean + value_se, colour = method), width = 0.2) + 
  facet_wrap(metric ~ ., scales = "free") + theme_bw()