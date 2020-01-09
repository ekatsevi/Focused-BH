FDP = matrix(0, num_eval_filters, num_methods)
power = matrix(0, num_eval_filters, num_methods)
nulls_retained = matrix(0, num_eval_filters, num_methods)
non_nulls_retained = matrix(0, num_eval_filters, num_methods)
rejections_retained = matrix(0, num_eval_filters, num_methods)
inflation_factor = matrix(0, num_eval_filters, num_methods)
avg_discovery_depth = matrix(0, num_eval_filters, num_methods)

for(method in 1:num_methods){
  R = rejections[method,]
  for(filter_idx in 1:num_eval_filters){
    eval_filter_type = eval_filter_types[filter_idx]
    if(eval_filter_type == "trivial"){
      filter_weights = R
    } else{
      filter_weights = get_filter_weights(R, G, eval_filter_type, base_weights)
    }
    V_filter = sum(filter_weights[nulls_node])
    S_filter = sum(filter_weights[!nulls_node])
    R_filter = sum(filter_weights)
    V_0 = sum(R[nulls_node])
    S_0 = sum(R[!nulls_node])
    R_0 = sum(R)
    
    FDP[filter_idx, method] = V_filter/R_filter
    power[filter_idx, method] = S_filter/max_true_rejections[filter_idx]
    nulls_retained[filter_idx, method] = V_filter/V_0
    non_nulls_retained[filter_idx, method] = S_filter/S_0
    rejections_retained[filter_idx, method] = R_filter/R_0
    avg_discovery_depth[filter_idx, method] = sum((G$depths)*filter_weights)/sum(filter_weights)
        
    if(R_filter == 0 & R_0 == 0){
      inflation_factor[filter_idx, method] = 1
    }
    if(R_filter == 0 & R_0 > 0){
      inflation_factor[filter_idx, method] = 0
    }
    if(R_filter > 0 & V_0 == 0){
      inflation_factor[filter_idx, method] = 1
    } 
    if(R_filter > 0 & V_0 > 0){
      inflation_factor[filter_idx, method] = (V_filter/R_filter)/(V_0/R_0)
    }
  }
}
FDP[is.na(FDP)] = 0
#nulls_retained[is.na(nulls_retained)] = 1
#non_nulls_retained[is.na(non_nulls_retained)] = 1
#rejections_retained[is.na(rejections_retained)] = 1

results = c()
results$FDP = FDP
results$power = power
results$nulls_retained = nulls_retained
results$non_nulls_retained = non_nulls_retained
results$rejections_retained = rejections_retained
results$inflation_factor = inflation_factor
results$avg_discovery_depth = avg_discovery_depth
