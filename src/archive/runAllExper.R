rm(list = ls()) # clean environment

params_version = 12
methods_version = 16
repeats = 500
#output_basename = "../output"
#precomp_basename = "../precomp"
filter_version = 3
eval_filter_type = "hard_outer_nodes"
eval_filter_types = c("trivial", "hard_outer_nodes")

B = 500 # number of resamples

# clear = TRUE
# clear_precomp_pvals = clear
# clear_precomp_gamma = clear
# clear_rejections = clear
# clear_evaluation = clear
clear_precomp_pvals = TRUE
clear_precomp_gamma = TRUE
clear_rejections = TRUE
clear_evaluation = TRUE


source("set_paths.R")
source("source_dependencies.R")

source("makeFilter.R")
source("makeMethods.R")
source("makeParams.R")
base_weights = make_base_weights(G, base_weights_type)
n = G$n
m = G$m

num_params = length(params)
num_eval_filters = length(eval_filter_types)

cat(sprintf("Running precomputation for filter %d, params version %d, methods version %d...\n", 
            filter_version, params_version, methods_version))

# create output directories if necessary
rejections_dir = sprintf("%s/rejections/p%d_m%d_f%d/", output_dir, params_version, methods_version, filter_version)
evaluation_dir = sprintf("%s/evaluation/p%d_m%d_f%d/", output_dir, params_version, methods_version, filter_version)
if(!file_test("-d", rejections_dir)){
  dir.create(rejections_dir)
}
if(!file_test("-d", evaluation_dir)){
  dir.create(evaluation_dir)
}

#Rprof(tmp <- tempfile(), line.profiling = TRUE)
if(precomp_dir != ""){
  source("precomp_local.R")
}

cat(sprintf("Running experiment for filter %d, params version %d, methods version %d...\n", 
            filter_version, params_version, methods_version))

for(param_idx in 1:num_params){
  cat(sprintf("Working on param_idx = %d...\n", param_idx))
  source("parse_params.R")
  
  # parameter-specific setup, if no precomputation has been done
  if(precomp_dir == ""){
    source("problem_setup_param.R") 
  } else{
    precomp_pvals_oracle_filename = sprintf("%s/pvals_p%d_%s.Rda", 
                                            precomp_dir, params_version, "oracle")
    precomp_gamma_perm_filename = sprintf("%s/gamma_p%d_f%d_%s.Rda", precomp_dir,
                                          params_version, filter_version, "permutation")
    precomp_gamma_oracle_filename = sprintf("%s/gamma_p%d_f%d_%s.Rda", precomp_dir,
                                            params_version, filter_version, "oracle")
    load(precomp_pvals_oracle_filename)
    load(precomp_gamma_perm_filename)
    load(precomp_gamma_oracle_filename)
  }
  
  for(rep in 1:repeats){
    rejections_filename = sprintf("%s/rejections/p%d_m%d_f%d/pi%d_r%d.Rda", output_dir, params_version, methods_version, 
                       filter_version, param_idx, rep)
    evaluation_filename = sprintf("%s/evaluation/p%d_m%d_f%d/pi%d_r%d.Rda", output_dir, params_version, methods_version, 
                                  filter_version, param_idx, rep)
    if(!file.exists(rejections_filename) | clear_rejections){
      if(precomp_dir == ""){
        source("problem_setup_rep.R")
      } else{
        P = P_oracle[rep,,param_idx]
      }
      source("run_methods.R")
      if(output_dir != ""){
        save(file = rejections_filename, rejections)
      }
    }
    
    if(!file.exists(evaluation_filename) | clear_evaluation){
      load(rejections_filename)
      source("evaluate_methods.R")
      if(output_dir != ""){
        save(file = evaluation_filename, results)
      }
    }
  }
}
#Rprof()
#summaryRprof(tmp, lines = "show")
#plot_experiment_4(params_version, methods_version, filter_version, eval_filter_types)  

source("plot_tree_experiment_job_talk.R")
