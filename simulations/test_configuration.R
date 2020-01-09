experiment_name = "PheWAS_Fisher"
input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
input_mode = "experiment"
source(input_filename, local = TRUE)

sum(nonnull_nodes)
mean(nonnull_nodes)

beta = numeric(num_items)
beta[nonnull_items] = 5
q = 0.1
filter_name = "outer_nodes"
FDP_vals = matrix(0,2,50, dimnames = list(c("before", "after"), NULL))
for(rep in 1:50){
  item_stats = beta + rnorm(num_items)
  P_item = pnorm(item_stats, lower.tail = FALSE)
  # P = sapply(ids, function(id)(min(sort(P_item[adj_matrix[,id]])*num_annotations[id]/(1:num_annotations[id]))))
  P = sapply(1:m, function(node)(pchisq(q = -2*sum(log(P_item[adj_matrix[,node]])), 
                                        df = 2*items_per_node[node], 
                                        lower.tail = FALSE)))
  R = p.adjust(P, "fdr") <= q # /mean(!nonnull_nodes)
  # df_output = REVIGO(P, R)
  # df_output = df_output %>% mutate(null = !nonnull_nodes[term_ID])
  R_F = run_filter(P, R, G, filter_name)
  print(sum(R_F[nonnull_nodes]))
  # FDP_vals[rep] = sum(R_F[!nonnull_nodes])/max(1,sum(R_F))
  # cat(sprintf("The FDP of BH is %0.2f for rep %d.\n", FDP_vals[rep], rep))
  # cat(sprintf("BH has %d/%d false positives before filtering (FDP = %0.2f), and %d/%d after filtering (FDP = %0.2f).\n",
  #             sum(R[!nonnull_nodes]), sum(R), sum(R[!nonnull_nodes])/sum(R),
  #             sum(R_F[!nonnull_nodes]), sum(R_F), sum(R_F[!nonnull_nodes])/sum(R_F)))
  # FDP_vals["before",rep] = sum(R[!nonnull_nodes])/sum(R)
  # FDP_vals["after",rep] = sum(R_F[!nonnull_nodes])/sum(R_F)
}
rowMeans(FDP_vals)
# df = melt(FDP_vals)
# names(df) = c("method", "filter", "rep", "FDP")
# df = as_tibble(df)
# cat(sprintf("The FDR of BH is %0.2f.\n", mean(FDP_vals)))

# df2 = tibble(before = FDP_vals["REVIGO", "before",],
#              after_REVIGO = FDP_vals["REVIGO", "after",],
#              after_REVIGO2 = FDP_vals["REVIGO2", "after",])
# df2 %>% summarise_all(mean)