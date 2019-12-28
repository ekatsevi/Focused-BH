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
  null_V_hats[[b]] = V_hat %>% filter(pvalue <= q)
}
null_V_hats = do.call("rbind", null_V_hats)
# null_V_hats = rbind(null_V_hats, tibble(pvalue = rep(0,B), V_oracle = rep(0,B), 
#                                         V_permutation = rep(0, B), b = 1:B))
# compute V_hat_permutation
V_hat_permutation = function(t){
  return(null_V_hats %>% filter(pvalue <= t) %>% group_by(b) %>% 
           summarise(V = max(V_permutation)) %>% summarise(mean(V)) %>% pull())
}


t = seq(0, 0.02, length.out = 1000)
df = tibble(t, V_hat_perm = sapply(t, V_hat_permutation), V_hat_original = m*t)
df %>% gather(estimate, V_hat, -t) %>%
  ggplot(aes(x = t, y = V_hat, group = estimate, colour = estimate)) +
  geom_line() + theme_bw()

