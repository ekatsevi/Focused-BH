run_one_experiment = function(experiment_name, experiment_index, base_dir){
  # set up all the simulation parameters
  
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  source(input_filename, local = TRUE)
  
  cat(sprintf("Starting simulation with signal strength = %0.1f...\n", signal_strength))
  
  # set up arrays to hold results
  num_methods = length(methods)
  
  # FDR and power
  metrics = array(0, 
                  dim = c(num_methods, 2, reps), 
                  dimnames = list(methods, c("power", "fdp"), NULL))

  # create set of non-nulls genes
  delta = 0
  beta = numeric(num_genes)
  beta[nonnull_genes] = signal_strength
  # logit_probabilities = base_log_odds + beta
  # gene_probabilities = plogis(logit_probabilities)
  
  # deduce non-null GO terms
  # genes_per_term = as.numeric(matrix(gene_probabilities,1,num_genes) %*% adj_matrix)
  # genes_per_term_null = mean(gene_probabilities)*colSums(adj_matrix)
  # nonnull_terms = genes_per_term > genes_per_term_null
  # means_per_term = as.numeric(matrix(beta,1,num_genes) %*% adj_matrix)/colSums(adj_matrix)
  # nonnull_terms = means_per_term > delta
  # mean(nonnull_terms)
  

  # mean(nonnull_terms)
  
  num_annotations = sapply(sets, length)
  
  # for permutation method
  gamma = function(t)(null_V_hats %>% filter(pvalue <= t) %>% summarise(mean(V_permutation)) %>% pull())
  gamma_oracle = function(t)(null_V_hats %>% filter(pvalue <= t) %>% summarise(mean(V_oracle)) %>% pull())
  
  for(rep in 1:reps){
    cat(sprintf("Working on repetition %d out of %d...\n", rep, reps))
    # problem setup

    DE_stats = beta + rnorm(num_genes)
    P_gene = pnorm(DE_stats, lower.tail = FALSE)
    P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])), 
                                        df = 2*num_annotations[id], 
                                        lower.tail = FALSE)))
    # P = numeric(m)
    # names(P) = ids
    # for(id in ids){
      # P[id] = fisher.test(DE_genes, adj_matrix[,id], alternative = "greater")$p.value
      # P[id] = t.test(DE_stats[adj_matrix[,id]], DE_stats[!adj_matrix[,id]], alternative = "greater")$p.value
    # }
    
    selected = matrix(FALSE, num_methods, m, dimnames = list(methods, ids))
    
    # BH
    selected["BH",] = p.adjust(P, "fdr") <= q/mean(!nonnull_terms)

    # Structured Holm
    if(min(P) >= q/m){
      selected["Structured_Holm",] = rep(FALSE, m)
    } else{
      selected["Structured_Holm",] = structuredHolm(G_SH, alpha_max = q, pvalues=P)@rejected
    }
        
    # Focused BH (original)
    selected["Focused_BH_original",] = FocusedBH(filter_name, P, G, q)
    
    # Focused BH (permutation)
    selected["Focused_BH_permutation",] = FocusedBH(filter_name, P, G, q, gamma)
    
    # Filter the results
    for(method in methods){
      R = selected[method,]
      R_F = run_filter(P, R, G, filter_name)
      # metrics["knockoffs", "power", rep] = sum(selected[nonnull_terms])/length(nonnull_terms)
      metrics[method, "power", rep] = sum(R_F[nonnull_terms])
      metrics[method, "fdp", rep] = sum(R_F[!nonnull_terms])/max(1,sum(R_F))
      cat(sprintf("Power of %s is %0.2f.\n", method, metrics[method,"power",rep]))
      cat(sprintf("FDP of %s is %0.2f.\n", method, metrics[method,"fdp",rep]))
    }
    
    # verbose option useful for interactive runs
    for(method in methods){
      cat(sprintf("Power of %s is %0.2f.\n", method, metrics[method,"power",rep]))
      cat(sprintf("FDP of %s is %0.2f.\n", method, metrics[method,"fdp",rep]))
    }
  }
  df = melt(metrics)
  names(df) = c("method", "metric", "rep", "value")
  df$signal_strength = signal_strength
  df = as_tibble(df)
  write_tsv(df, path = sprintf("%s/results/%s/metrics_%s_%d.tsv", base_dir, experiment_name, experiment_name, experiment_index))
}

run_one_precomputation = function(experiment_name, b, base_dir){
  set.seed(1234+b) # for reproducibility
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  source(input_filename, local = TRUE)
  
  num_annotations = sapply(sets, length)
  
  cat(sprintf("Generating null p-values...\n"))
  DE_stats = rnorm(num_genes)
  P_gene = pnorm(DE_stats, lower.tail = FALSE)
  P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])), 
                                      df = 2*num_annotations[id], 
                                      lower.tail = FALSE)))
  ord = order(P)
  
  k_max = sum(P <= t_max)
  V_hat = matrix(0, k_max, 3, dimnames = list(NULL, c("pvalue", "V_oracle", "V_permutation")))
  
  
  for(k in k_max:1){
    cat(sprintf("Considering rejection set of size %d...\n", k))
    R = logical(m)
    R[ord[1:k]] = TRUE
    names(R) = ids
    R_F = run_filter(P, R, G, filter_name)
    V_hat[k, "pvalue"] = P[ord[k]]
    V_hat[k, "V_oracle"] = sum(R_F & !nonnull_terms)
    V_hat[k, "V_permutation"] = sum(R_F)
  }
  write_tsv(as_tibble(V_hat), 
            path = sprintf("%s/precomp/%s/precomp_%s_%d.tsv", 
                           base_dir, experiment_name, experiment_name, b))
}