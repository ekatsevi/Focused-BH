# precompute the permutation estimate V_hat
run_one_precomputation = function(experiment_name, experiment_index, base_dir){
  set.seed(1234+experiment_index) # for reproducibility; should be different for different experiment_index
  
  # set up all the simulation parameters
  input_filename = sprintf("simulations/input_files/input_file_%s.R", experiment_name)
  input_mode = "precomputation"
  source(input_filename, local = TRUE)
  
  # generate the null p-values
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
  
  # matrix to store V_hat 
  V_hat = matrix(0, k_max+1, 2, dimnames = list(NULL, c("pvalue", "V_permutation")))
  
  # treat outer nodes filter separately because a faster algorithm is available
  if(filter_name == "outer_nodes"){
    V_hat[,"pvalue"] = c(0,sort(P))
    V_hat[,"V_permutation"] = rev(get_num_outer_nodes(P, G))
  } else{
    for(k in k_max:1){
      cat(sprintf("Considering rejection set of size %d...\n", k))
      R = logical(m)
      R[ord[1:k]] = TRUE
      R_F = run_filter(P, R, G, filter_name)
      V_hat[k+1, "pvalue"] = P[ord[k]]      
      V_hat[k+1, "V_permutation"] = sum(R_F)
    }
    V_hat[1, "pvalue"] = 0
    V_hat[1, "V_permutation"] = 0
  }
  
  # write to file
  output_filename = sprintf("%s/precomp/%s/precomp_%s_%d.tsv", 
                           base_dir, experiment_name, experiment_name, b) 
  write_tsv(as_tibble(V_hat), path = output_filename)
  
  # add provenance information to output file
  call_info = sprintf("#\n#\n### FUNCTION CALL: ###\n# run_one_precomputation(experiment_name = \"%s\", experiment_index = %d, base_dir = \"%s\")", 
                      experiment_name, experiment_index, base_dir)
  write(call_info, output_filename, append = TRUE)
  add_git_info(output_filename)
  add_R_session_info(output_filename)
}