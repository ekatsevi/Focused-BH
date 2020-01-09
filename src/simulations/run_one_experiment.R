# run one experiment
run_one_experiment = function(experiment_name, experiment_index, base_dir){
  # set up all the simulation parameters
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  input_mode = "experiment"
  source(input_filename, local = TRUE)
  
  # set up arrays to hold results (FDP and power)
  num_methods = length(methods)
  num_parameters = nrow(parameters)
  metrics = array(0, 
                  dim = c(num_methods, 2, num_parameters), 
                  dimnames = list(methods, c("power", "fdp"), NULL))
  
  # define V_hat for Focused BH (permutation) based on precomputed null values
  if("Focused_BH_permutation" %in% methods){
    # load precomputation results
    precomp_dir = sprintf("%s/precomp/%s", base_dir, experiment_name)
    stopifnot(dir.exists(precomp_dir))
    precomp_files = list.files(precomp_dir)
    B = length(precomp_files)
    stopifnot(B > 0)
    null_V_hats = vector("list", B)
    for(b in 1:B){
      precomp_filename = sprintf("%s/precomp/%s/%s",base_dir, experiment_name, precomp_files[b])
      V_hat = read_tsv(precomp_filename, comment = "#")
      V_hat$b = b
      null_V_hats[[b]] = V_hat %>% filter(pvalue <= q)
    }
    null_V_hats = do.call("rbind", null_V_hats)
    
    # compute V_hat_permutation
    V_hat_permutation = function(t){
      return(null_V_hats %>% filter(pvalue <= t) %>% group_by(b) %>% 
               summarise(V = max(V_permutation)) %>% summarise(mean(V)) %>% pull())
    }
  }
  
  # precomputation for Yekutieli
  if("Yekutieli" %in% methods){
    depths = get_depths(G$Pa)
    q_Yekutieli = q/(2*max(depths)) # corrected FDR target level based on 
                                    # Proposition 1 of Yekutieli (2008)
    
  }
  
  # run all methods for each setting of parameters 
  for(index in 1:num_parameters){
    rep = parameters$rep[index] 
    signal_strength = parameters$signal_strength[index]
    
    set.seed(1234+1000*experiment_index + rep) # for reproducibility; should vary by rep
    
    # create effect sizes based on signal strength and non-null items
    beta = numeric(num_items)
    beta[nonnull_items] = signal_strength
    
    cat(sprintf("Working on parameter index %d out of %d, corresponding to rep %d for signal strength %0.1f...\n", 
                index, num_parameters, rep, signal_strength))
    
    ### Problem setup ###
    # compute item-level p-values
    item_stats = beta + rnorm(num_items)
    P_item = pnorm(item_stats, lower.tail = FALSE)
    
    # compute node-level p-values
    if(global_test == "Fisher"){
      P = sapply(1:m, function(node)(pchisq(q = -2*sum(log(P_item[adj_matrix[,node]])), 
                                            df = 2*items_per_node[node], 
                                            lower.tail = FALSE)))
    }
    if(global_test == "Simes"){
      P = sapply(1:m, function(node)(min(sort(P_item[adj_matrix[,node]])*items_per_node[node]/(1:items_per_node[node]))))
    }
    
    ### Run all methods ###
    
    # to store rejections
    rejections = matrix(FALSE, num_methods, m, dimnames = list(methods, 1:m))
    
    # BH
    if("BH" %in% methods){
      rejections["BH",] = p.adjust(P, "fdr") <= q
    }
    
    # Storey-BH
    if("Storey_BH" %in% methods){
      rejections["Storey_BH",] = p.adjust(P, "fdr") <= q/mean(!nonnull_terms)
    }
    
    # Structured Holm
    if("Structured_Holm" %in% methods){
      if(min(P) >= q/m){
        rejections["Structured_Holm",] = rep(FALSE, m)
      } else{
        rejections["Structured_Holm",] = structuredHolm(G_SH, alpha_max = q, pvalues=P)@rejected
      }
    }
    
    if("Yekutieli" %in% methods){
      rejections["Yekutieli",] = Yekutieli(P, G_Yekutieli, q_Yekutieli)
    }
    
    if(exists("node_heights")){
      for(height in 1:max(node_heights)){
        leaf_bh_method = sprintf("Leaf_BH_%d", height)
        if(leaf_bh_method %in% methods){
          nodes_at_height = node_heights == height
          rejections[leaf_bh_method,nodes_at_height] = p.adjust(P[nodes_at_height], "fdr") <= q
        }
      }
    }
    
    fbh_methods = intersect(methods, c("Focused_BH_original", 
                                       "Focused_BH_permutation", 
                                       "Focused_BH_oracle"))
    
    if(length(fbh_methods) > 0){
      V_hats = vector("list", length(fbh_methods))
      names(V_hats) = fbh_methods
      
      # Focused BH (original)
      if("Focused_BH_original" %in% methods){
        V_hats[["Focused_BH_original"]] = function(t)(m*t)
      }
      
      # Focused BH (permutation)
      if("Focused_BH_permutation" %in% methods){
        V_hats[["Focused_BH_permutation"]] = V_hat_permutation
      }
      
      rejections[fbh_methods,] = FocusedBH(filter_name, P, G, q, V_hats)
    }
    
    # Filter the results
    for(method in methods){
      R = rejections[method,]
      if(sum(R) > k_max){
        R[which(R)[order(P[R])[(k_max+1):sum(R)]]] = FALSE
      }
      R_F = run_filter(P, R, G, filter_name)
      metrics[method, "power", index] = sum(R_F[nonnull_nodes])
      metrics[method, "fdp", index] = sum(R_F[!nonnull_nodes])/max(1,sum(R_F))
      cat(sprintf("Power of %s is %0.2f.\n", method, metrics[method,"power",index]))
      cat(sprintf("FDP of %s is %0.2f.\n", method, metrics[method,"fdp",index]))
    }
  }
  
  # massage results into data frame and write to file
  df = melt(metrics)
  names(df) = c("method", "metric", "index", "value")
  df = as_tibble(df)
  df = df %>% mutate(signal_strength = parameters$signal_strength[index], 
                     rep = parameters$rep[index]) %>%
    dplyr::select(-index)
  output_filename = sprintf("%s/results/%s/metrics_%s_%d.tsv", 
                            base_dir, experiment_name, experiment_name, experiment_index)
  write_tsv(df, path = output_filename)
  
  # add provenance information to output file
  call_info = sprintf("#\n#\n### FUNCTION CALL: ###\n# run_one_experiment(experiment_name = \"%s\", experiment_index = %d, base_dir = \"%s\")", 
                      experiment_name, experiment_index, base_dir)
  write(call_info, output_filename, append = TRUE)
  add_git_info(output_filename)
  add_R_session_info(output_filename)
}
