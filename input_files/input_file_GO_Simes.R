suppressPackageStartupMessages(library(tidyverse))

### set top-level parameters
reps = 25                         # Number of outer-loop repetitions
# reps = 1
methods = c("BH",                 # methods
            "Focused_BH_original",   
            "Focused_BH_permutation",
            # "Focused_BH_oracle",
            "Structured_Holm")
q = 0.1                           # FDR control level
signal_strength_vals = seq(2,6,by=0.5)
# signal_strength_vals = 6
global_test = "Simes"
k_max = 350   # maximum p-value threshold considered by permutation approach
B = 100       # number of repetitions
filter_name = "REVIGO"
parameters = tibble(signal_strength_vals)
names(parameters) = "signal_strength"

### set input_mode (experiment, precomputation, num_experiments, num_precomputations)
if(exists("input_mode")){
  stopifnot(input_mode %in% c("experiment", "precomputation"))
  if(input_mode == "experiment"){
    stopifnot(exists("experiment_index"))
    signal_strength = parameters$signal_strength[experiment_index]
  }
  if(input_mode == "precomputation"){
    stopifnot(exists("b"))
  }
} else{
  args = commandArgs(trailingOnly = TRUE)
  num_args = length(args)
  stopifnot(num_args == 1)
  input_mode = args[1]
  stopifnot(input_mode %in% c("num_experiments", "num_precomputations"))
  if(input_mode == "num_experiments"){
    cat(nrow(parameters))
  }
  if(input_mode == "num_precomputations"){
    cat(B)
  }
}

### define the graph
if(input_mode %in% c("precomputation", "experiment")){
  # read in GO data
  cat(sprintf("Reading in GO data...\n"))
  reduced_graph_100_file = sprintf("%s/data/processed/reduced_graph_100.Rda", base_dir)
  load(reduced_graph_100_file)
  # reduced_graph_file = sprintf("%s/data/processed/reduced_graph.Rda", base_dir)
  # load(reduced_graph_file)
  
  gene_sets = G$gene_sets
  ids = names(G$C)
  m = G$m
  num_annotations = sapply(gene_sets, length)
  
  index = c()
  for(i in 1:m){
    index[[ids[i]]] = i
  }
  
  parents_indexed = sapply(G$Pa, function(parents_list)(unname(index[parents_list])))
  children_indexed = sapply(G$C, function(children_list)(unname(index[children_list])))

  G = c()
  G$m = m
  G$Pa = parents_indexed
  G$C = children_indexed
  
  if("Structured_Holm" %in% methods & input_mode == "experiment"){
    cat(sprintf("Creating DAG for Structured Holm...\n"))
    G_SH = new("DAGstructure", parents = parents_indexed, 
               children = children_indexed, sets = gene_sets, twoway = FALSE)
  }
  
  genes = unique(unlist(gene_sets))
  num_genes = length(genes)
  adj_matrix = matrix(FALSE, num_genes, m)
  rownames(adj_matrix) = genes
  colnames(adj_matrix) = ids
  for(id in ids){
    adj_matrix[gene_sets[[id]],id] = TRUE
  }
  
  # anchor_terms = withSeed(sample(which(num_annotations == 5), 3), 1)
  # nonnull_genes = unique(unlist(lapply(gene_sets[anchor_terms], function(set)(set))))
  # nonnull_terms = colSums(adj_matrix[nonnull_genes,]) > 0

  anchor_terms = withSeed(sample(which(num_annotations == 5), 2), 1)
  nonnull_genes = unique(unlist(lapply(gene_sets[anchor_terms], function(set)(set))))
  nonnull_terms = colSums(adj_matrix[nonnull_genes,]) > 0
  
  # sum(nonnull_terms)
  # mean(nonnull_terms)
  
  # beta = numeric(num_genes)
  # names(beta) = genes
  # beta[nonnull_genes] = Inf
  # q = 0.1
  # k_max = 350
  # filter_name = "REVIGO"
  # FDP_vals = array(0,c(2,2,50), list(c("REVIGO", "REVIGO2"), c("before", "after"), NULL))
  # for(rep in 1:50){
  #   DE_stats = beta + rnorm(num_genes)
  #   P_gene = pnorm(DE_stats, lower.tail = FALSE)
  #   P = sapply(ids, function(id)(min(sort(P_gene[adj_matrix[,id]])*num_annotations[id]/(1:num_annotations[id]))))
  #   # P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])),
  #   #                                     df = 2*num_annotations[id],
  #   #                                     lower.tail = FALSE)))
  #   R = p.adjust(P, "fdr") <= q/mean(!nonnull_terms)
  #   if(sum(R) > k_max){
  #     R[which(R)[order(P[R])[(k_max+1):sum(R)]]] = FALSE
  #   }
  #   # df_output = REVIGO(P, R)
  #   # df_output = df_output %>% mutate(null = !nonnull_terms[term_ID])
  #   for(filter_name in c("REVIGO")){
  #     # cat(sprintf("Working on filter %s...\n", filter_name))
  #     R_F = run_filter(P, R, G, filter_name)
  #     # FDP_vals[rep] = sum(R_F[!nonnull_terms])/max(1,sum(R_F))
  #     # cat(sprintf("The FDP of BH is %0.2f for rep %d.\n", FDP_vals[rep], rep))
  #     cat(sprintf("BH has %d/%d false positives before filtering (FDP = %0.2f), and %d/%d after filtering (FDP = %0.2f).\n", 
  #                 sum(R[!nonnull_terms]), sum(R), sum(R[!nonnull_terms])/sum(R),
  #                 sum(R_F[!nonnull_terms]), sum(R_F), sum(R_F[!nonnull_terms])/sum(R_F)))
  #     FDP_vals[filter_name, "before",rep] = sum(R[!nonnull_terms])/sum(R)
  #     FDP_vals[filter_name, "after",rep] = sum(R_F[!nonnull_terms])/sum(R_F)
  #   }
  # }
  # df = melt(FDP_vals)
  # names(df) = c("method", "filter", "rep", "FDP")
  # df = as_tibble(df)
  # # cat(sprintf("The FDR of BH is %0.2f.\n", mean(FDP_vals)))
  # 
  # df2 = tibble(before = FDP_vals["REVIGO", "before",], 
  #              after_REVIGO = FDP_vals["REVIGO", "after",],
  #              after_REVIGO2 = FDP_vals["REVIGO2", "after",])
  # df2 %>% summarise_all(mean)
  # 
  # df2 %>% gather(filter_name, after_FDP, -before) %>% dplyr::rename(before_FDP = before) %>%
  #   ggplot(aes(x = before_FDP, y = after_FDP, group = filter_name, colour = filter_name)) + 
  #   geom_point() + theme_bw() + geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed")
  # 
  # #######################################
  # beta = numeric(num_genes)
  # names(beta) = genes
  # beta[nonnull_genes] = 5
  # q = 0.1
  # reps = 25
  # FDP_vals = numeric(reps)
  # for(rep in 1:reps){
  #   DE_stats = beta + rnorm(num_genes)
  #   P_gene = pnorm(DE_stats, lower.tail = FALSE)
  #   ord = order(P_gene)
  #   P_gene_sorted = P_gene[ord]
  #   adj_matrix_sorted = adj_matrix[ord,]
  #   P = sapply(ids, function(id)(min(sort(P_gene
  #                                         [adj_matrix[,id]])*num_annotations[id]/(1:num_annotations[id]))))
  #   # P = sapply(ids, function(id)(min(P_gene_sorted[adj_matrix_sorted[,id]]*num_annotations[id]/(1:num_annotations[id]))))
  #   # P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])),
  #   #                                     df = 2*num_annotations[id],
  #   #                                     lower.tail = FALSE)))
  #   R = p.adjust(P, "fdr") <= q/mean(!nonnull_terms)
  #   print(sum(R[nonnull_terms]))
  #   # FDP_vals[rep] = sum(R[!nonnull_terms])/sum(R)
  #   # print(FDP_vals[rep])
  # }

  
  # sum(nonnull_terms)
  # mean(nonnull_terms)
  # P = numeric(m)
  # P[!nonnull_terms] = 1
  # R_F = run_filter(P, nonnull_terms, G, "REVIGO")
  # hist(colSums(adj_matrix[nonnull_genes,nonnull_terms]), breaks = 10)
  # hist(colSums(adj_matrix[nonnull_genes,nonnull_terms])/colSums(adj_matrix[,nonnull_terms]), breaks = 20)
}