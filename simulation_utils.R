run_one_experiment = function(experiment_name, experiment_index, base_dir){
  set.seed(1234) # for reproducibility
  
  # set up all the simulation parameters
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  mode = "experiment"
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

  num_annotations = sapply(sets, length)
  
  # for permutation method
  if("Focused_BH_permutation" %in% methods){
    # load precomputation results
    precomp_dir = sprintf("%s/precomp/%s", base_dir, experiment_name)
    stopifnot(dir.exists(precomp_dir))
    precomp_files = list.files(precomp_dir)
    B = length(precomp_files)
    stopifnot(B > 0)
    null_V_hats = vector("list", length())
    for(b in 1:B){
        V_hat = read_tsv(precomp_files[b], col_types = "dii")
        V_hat$b = b
        null_V_hats[[b]] = V_hat
    }
    null_V_hats = do.call("rbind", null_V_hats)
    # compute V_hat_permutation
    V_hat_permutation = function(t){
      if(t <= t_max){
        return(null_V_hats %>% filter(pvalue <= t) %>% summarise(mean(V_permutation)) %>% pull())
      } else{
        return(m*t)
      }
    }
    # compute V_hat_oracle
    V_hat_oracle = function(t){
      if(t <= t_max){
        return(null_V_hats %>% filter(pvalue <= t) %>% summarise(mean(V_oracle)) %>% pull())
      } else{
        return(m*t)
      }
    }
  }
  
  for(rep in 1:reps){
    cat(sprintf("Working on repetition %d out of %d...\n", rep, reps))
    
    ### Problem setup ###
    # compute gene-level p-values
    DE_stats = beta + rnorm(num_genes)
    P_gene = pnorm(DE_stats, lower.tail = FALSE)
    # compute group-level p-values
    if(global_test == "Fisher"){
      P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])), 
                                          df = 2*num_annotations[id], 
                                          lower.tail = FALSE)))
    }
    if(global_test == "Simes"){
      P = sapply(ids, function(id)(min(sort(P_gene[adj_matrix[,id]])*num_annotations/(1:num_annotations))))
    }
    
    ### Run all methods ###
    
    # to store rejections
    rejections = matrix(FALSE, num_methods, m, dimnames = list(methods, ids))
    
    # BH
    if("BH" %in% methods){
      rejections["BH",] = p.adjust(P, "fdr") <= q/mean(!nonnull_terms)
    }

    # Structured Holm
    if("Structured_Holm" %in% methods){
      if(min(P) >= q/m){
        rejections["Structured_Holm",] = rep(FALSE, m)
      } else{
        rejections["Structured_Holm",] = structuredHolm(G_SH, alpha_max = q, pvalues=P)@rejected
      }
    }
        
    # Focused BH (original)
    if("Focused_BH_original" %in% methods){
      rejections["Focused_BH_original",] = FocusedBH(filter_name, P, G, q)
    }
    
    # Focused BH (permutation)
    if("Focused_BH_permutation" %in% methods){
      rejections["Focused_BH_permutation",] = FocusedBH(filter_name, P, G, q, V_hat_permutation)      
    }
    
    # Focused BH (oracle)
    if("Focused_BH_oracle" %in% methods){
      rejections["Focused_BH_oracle",] = FocusedBH(filter_name, P, G, q, V_hat_oracle)      
    }
    
    # Filter the results
    for(method in methods){
      R = rejections[method,]
      R_F = run_filter(P, R, G, filter_name)
      # metrics["knockoffs", "power", rep] = sum(rejections[nonnull_terms])/length(nonnull_terms)
      metrics[method, "power", rep] = sum(R_F[nonnull_terms])
      metrics[method, "fdp", rep] = sum(R_F[!nonnull_terms])/max(1,sum(R_F))
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
  set.seed(1234+b) # for reproducibility; should be different for different b
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  mode = "precomputation"
  source(input_filename, local = TRUE)
  
  cat(sprintf("Generating null p-values...\n"))
  DE_stats = rnorm(num_genes)
  P_gene = pnorm(DE_stats, lower.tail = FALSE)
  if(global_test == "Fisher"){
    P = sapply(ids, function(id)(pchisq(q = -2*sum(log(P_gene[adj_matrix[,id]])), 
                                        df = 2*num_annotations[id], 
                                        lower.tail = FALSE)))
  }
  if(global_test == "Simes"){
    P = sapply(ids, function(id)(min(sort(P_gene[adj_matrix[,id]])*num_annotations[id]/(1:num_annotations[id]))))
  }
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