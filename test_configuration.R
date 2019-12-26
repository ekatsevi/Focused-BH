experiment_name = "GO_Simes"
input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
input_mode = "experiment"
source(input_filename, local = TRUE)

sum(nonnull_terms)
mean(nonnull_terms)

beta = numeric(num_genes)
names(beta) = genes
beta[nonnull_genes] = Inf
q = 0.1
k_max = 350
filter_name = "REVIGO"
FDP_vals = array(0,c(2,2,50), list(c("REVIGO", "REVIGO2"), c("before", "after"), NULL))
for(rep in 1:50){
  DE_stats = beta + rnorm(num_genes)
  P_gene = pnorm(DE_stats, lower.tail = FALSE)
  P = sapply(ids, function(id)(min(sort(P_gene[adj_matrix[,id]])*num_annotations[id]/(1:num_annotations[id]))))
  # P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])),
  #                                     df = 2*num_annotations[id],
  #                                     lower.tail = FALSE)))
  R = p.adjust(P, "fdr") <= q/mean(!nonnull_terms)
  if(sum(R) > k_max){
    R[which(R)[order(P[R])[(k_max+1):sum(R)]]] = FALSE
  }
  # df_output = REVIGO(P, R)
  # df_output = df_output %>% mutate(null = !nonnull_terms[term_ID])
  for(filter_name in c("REVIGO")){
    # cat(sprintf("Working on filter %s...\n", filter_name))
    R_F = run_filter(P, R, G, filter_name)
    # FDP_vals[rep] = sum(R_F[!nonnull_terms])/max(1,sum(R_F))
    # cat(sprintf("The FDP of BH is %0.2f for rep %d.\n", FDP_vals[rep], rep))
    cat(sprintf("BH has %d/%d false positives before filtering (FDP = %0.2f), and %d/%d after filtering (FDP = %0.2f).\n",
                sum(R[!nonnull_terms]), sum(R), sum(R[!nonnull_terms])/sum(R),
                sum(R_F[!nonnull_terms]), sum(R_F), sum(R_F[!nonnull_terms])/sum(R_F)))
    FDP_vals[filter_name, "before",rep] = sum(R[!nonnull_terms])/sum(R)
    FDP_vals[filter_name, "after",rep] = sum(R_F[!nonnull_terms])/sum(R_F)
  }
}
df = melt(FDP_vals)
names(df) = c("method", "filter", "rep", "FDP")
df = as_tibble(df)
# cat(sprintf("The FDR of BH is %0.2f.\n", mean(FDP_vals)))

df2 = tibble(before = FDP_vals["REVIGO", "before",],
             after_REVIGO = FDP_vals["REVIGO", "after",],
             after_REVIGO2 = FDP_vals["REVIGO2", "after",])
df2 %>% summarise_all(mean)