run_one_experiment = function(experiment_name, experiment_index, base_dir){
  
  # set up all the simulation parameters
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  input_mode = "experiment"
  source(input_filename, local = TRUE)

  # set up arrays to hold results
  num_methods = length(methods)
  num_parameters = nrow(parameters)
  # FDR and power
  metrics = array(0, 
                  dim = c(num_methods, 2, num_parameters), 
                  dimnames = list(methods, c("power", "fdp"), NULL))

  # for permutation method
  if("Focused_BH_permutation" %in% methods){
    # load precomputation results
    precomp_dir = sprintf("%s/precomp/%s", base_dir, experiment_name)
    stopifnot(dir.exists(precomp_dir))
    precomp_files = list.files(precomp_dir)
    B = length(precomp_files)
    stopifnot(B > 0)
    null_V_hats = vector("list", B)
    for(b in 1:B){
        V_hat = read_tsv(sprintf("%s/precomp/%s/%s",base_dir, experiment_name, precomp_files[b]), col_types = "dii")
        V_hat$b = b
        null_V_hats[[b]] = V_hat
    }
    null_V_hats = do.call("rbind", null_V_hats)
    null_V_hats = rbind(null_V_hats, tibble(pvalue = rep(0,B), V_oracle = rep(0,B), 
                                               V_permutation = rep(0, B), b = 1:B))
    # compute V_hat_permutation
    V_hat_permutation = function(t){
      return(null_V_hats %>% filter(pvalue <= t) %>% group_by(b) %>% 
               summarise(V = max(V_permutation)) %>% summarise(mean(V)) %>% pull())
    }
    # compute V_hat_oracle
    V_hat_oracle = function(t){
      return(null_V_hats %>% filter(pvalue <= t) %>% group_by(b) %>% 
               summarise(V = max(V_oracle)) %>% summarise(mean(V)) %>% pull())
    }
  }
  
  for(index in 1:num_parameters){
    rep = parameters$rep[index]
    signal_strength = parameters$signal_strength[index]

    set.seed(1234+rep) # for reproducibility; should vary by rep
    
    # create set of non-nulls genes
    beta = numeric(num_items)
    beta[nonnull_items] = signal_strength
    
    
    cat(sprintf("Working on parameter index %d out of %d, corresponding to rep %d for signal strength %0.1f...\n", 
                index, num_parameters, rep, signal_strength))
    
    ### Problem setup ###
    # compute gene-level p-values
    item_stats = beta + rnorm(num_items)
    P_item = pnorm(item_stats, lower.tail = FALSE)
    # compute group-level p-values
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
        
    fbh_methods = intersect(methods, c("Focused_BH_original", 
                                       "Focused_BH_permutation", 
                                       "Focused_BH_oracle"))
    
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
    
    # Focused BH (oracle)
    if("Focused_BH_oracle" %in% methods){
      V_hats[["Focused_BH_oracle"]] = V_hat_oracle
    }
    
    rejections[fbh_methods,] = FocusedBH(filter_name, P, G, q, V_hats)
    
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
  df = melt(metrics)
  names(df) = c("method", "metric", "index", "value")
  df = as_tibble(df)
  df = df %>% mutate(signal_strength = parameters$signal_strength[index], rep = parameters$rep[index]) %>%
    dplyr::select(-index)
  output_filename = sprintf("%s/results/%s/metrics_%s_%d.tsv", 
                            base_dir, experiment_name, experiment_name, experiment_index)
  write_tsv(df, path = output_filename)
  
  # add provenance information to output file
  add_git_info(output_filename)
  add_R_session_info(output_filename)
}

run_one_precomputation = function(experiment_name, b, base_dir){
  set.seed(1234+b) # for reproducibility; should be different for different b
  input_filename = sprintf("input_files/input_file_%s.R", experiment_name)
  input_mode = "precomputation"
  source(input_filename, local = TRUE)
  
  cat(sprintf("Generating null p-values...\n"))
  item_stats = rnorm(num_items)
  P_item = pnorm(item_stats, lower.tail = FALSE)
  if(global_test == "Fisher"){
    P = sapply(1:m, function(node)(pchisq(q = -2*sum(log(P_item[adj_matrix[,node]])), 
                                        df = 2*items_per_node[node], 
                                        lower.tail = FALSE)))
  }
  if(global_test == "Simes"){
    P = sapply(1:m, function(node)(min(sort(P_item[adj_matrix[,node]])*items_per_node[node]/(1:items_per_node[node]))))
  }
  ord = order(P)
  
  V_hat = matrix(0, k_max, 3, dimnames = list(NULL, c("pvalue", "V_oracle", "V_permutation")))
  
  for(k in k_max:1){
    cat(sprintf("Considering rejection set of size %d...\n", k))
    R = logical(m)
    R[ord[1:k]] = TRUE
    R_F = run_filter(P, R, G, filter_name)
    V_hat[k, "pvalue"] = P[ord[k]]
    V_hat[k, "V_oracle"] = sum(R_F & !nonnull_terms)
    V_hat[k, "V_permutation"] = sum(R_F)
  }
  write_tsv(as_tibble(V_hat), 
            path = sprintf("%s/precomp/%s/precomp_%s_%d.tsv", 
                           base_dir, experiment_name, experiment_name, b))
}