######################################################
# CREATE PLOT COMPARING PERMUTATION TO ORIGINAL V_HAT 
######################################################

# setup
machine = "local"
source("setup.R")

# extract permutation V_hat for GO/REVIGO experiment
experiment_name = "GO_Simes"
GO_reduced_graph_100_file = sprintf("%s/data/processed/GO/GO_reduced_graph_100.Rda", base_dir)
load(GO_reduced_graph_100_file)
m_GO = G$m
precomp_dir = sprintf("%s/precomp/%s", base_dir, experiment_name)
stopifnot(dir.exists(precomp_dir))
precomp_files = list.files(precomp_dir)
B = length(precomp_files)
stopifnot(B > 0)
null_V_hats = vector("list", B)
for(b in 1:B){
  V_hat = read_tsv(sprintf("%s/precomp/%s/%s", base_dir, experiment_name, precomp_files[b]), 
                   col_types = "dii", comment = "#")
  V_hat$b = b
  null_V_hats[[b]] = V_hat
}
null_V_hats_GO = do.call("rbind", null_V_hats)
V_hat_permutation_GO = function(t){
  return(null_V_hats_GO %>% filter(pvalue <= t) %>% group_by(b) %>% 
           summarise(V = max(V_permutation)) %>% summarise(mean(V)) %>% pull())
}

# extract permutation V_hat for ICD/outer nodes experiment
experiment_name = "PheWAS_Fisher"
graph_file = sprintf("%s/data/processed/biobank/ICD_graph.Rda", base_dir)
load(graph_file)
m_PheWAS = G$m
precomp_dir = sprintf("%s/precomp/%s", base_dir, experiment_name)
stopifnot(dir.exists(precomp_dir))
precomp_files = list.files(precomp_dir)
B = length(precomp_files)
stopifnot(B > 0)
null_V_hats = vector("list", B)
for(b in 1:B){
  V_hat = read_tsv(sprintf("%s/precomp/%s/%s", base_dir, experiment_name, precomp_files[b]), 
                   col_types = "di", comment = "#")
  V_hat$b = b
  null_V_hats[[b]] = V_hat %>% filter(pvalue <= 0.05)
}
null_V_hats_PheWAS = do.call("rbind", null_V_hats)
V_hat_permutation_PheWAS = function(t){
  return(null_V_hats_PheWAS %>% filter(pvalue <= t) %>% group_by(b) %>% 
           summarise(V = max(V_permutation)) %>% summarise(mean(V)) %>% pull())
}

# combine results together and plot
t = seq(0, 0.02, length.out = 1000)
df = tibble(t, 
            V_hat_GO = sapply(t, V_hat_permutation_GO)/m_GO, 
            V_hat_PheWAS = sapply(t, V_hat_permutation_PheWAS)/m_PheWAS, 
            V_hat_original = t)
p = df %>% gather(estimate, V_hat, -t) %>%
  mutate(estimate = factor(estimate, 
                           levels = c("V_hat_original", "V_hat_PheWAS", "V_hat_GO"),
                           labels = c("Original", "Permutation (ICD)", "Permutation (GO)"))) %>%
  ggplot(aes(x = t, y = V_hat, group = estimate, colour = estimate, linetype = estimate)) +
  scale_colour_manual(values = c("dodgerblue", "blue", "blue")) + 
  scale_linetype_manual(values = c("solid", "longdash", "dotdash")) + 
  geom_line() + theme_bw() + theme(legend.position = c(0.275, 0.8), legend.title = element_blank()) + 
  xlab("p-value threshold") + ylab("Estimated false discoveries per hypothesis") + coord_fixed()
plot(p)
plot_filename = sprintf("%s/permutation_estimates.pdf", figures_dir) 
ggsave(filename = plot_filename, plot = p, device = "pdf", width = 3.75, height = 3.75)