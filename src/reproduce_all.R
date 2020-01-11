################################################################
# REPRODUCE ALL NUMERICAL SIMULATIONS AND DATA ANALYSIS RESULTS
################################################################

# clear workspace
rm(list = ls())

# set-up, pre-processing
machine = "local"
source("setup.R")

############# RUN NUMERICAL SIMULATIONS ###################

# precomputation for GO_Simes, PheWAS_intermediate simulation
experiment_names = c("GO_Simes", "PheWAS_intermediate")
for(experiment_name in experiment_names){
  cat(sprintf("Running precomputation for experiment %s...\n", experiment_name))
  precomp_dir=sprintf("%s/precomp/%s", base_dir, experiment_name)
  if(!dir.exists(precomp_dir)){
    dir.create(precomp_dir)
  }
  input_mode = "num_precomputations"
  B = as.numeric(captureOutput(source(sprintf("simulations/input_files/input_file_%s.R", experiment_name))))
  for(b in 1:B){
    cat(sprintf("Computing permutation number %d out of %d...\n", b, B))
    run_one_precomputation(experiment_name, b, base_dir)  
  }
}

# numerical simulations
experiment_names = c("GO_Simes", "PheWAS_dispersed", "PheWAS_intermediate", 
                     "PheWAS_clustered", "PheWAS_intermediate_leaf")
for(experiment_name in experiment_names){
  cat(sprintf("Running experiment %s...\n", experiment_name))
  results_dir=sprintf("%s/results/%s", base_dir, experiment_name)
  if(!dir.exists(results_dir)){
    dir.create(results_dir)
  }
  input_mode = "num_experiments"
  num_experiments = as.numeric(captureOutput(source(sprintf("simulations/input_files/input_file_%s.R", experiment_name))))
  for(experiment_index in 1:num_experiments){
    cat(sprintf("Running experiment number %d out of %d...\n", experiment_index, num_experiments))
    run_one_experiment(experiment_name, experiment_index, base_dir)  
  }
}

############# RUN DATA ANALYSES ###################

# PheWAS data analysis
source("data_analysis/PheWAS_data_analysis.R")

# REVIGO data analysis
source("data_analysis/REVIGO_data_analysis.R")

############# PLOT RESULTS ###################

# plot figures for numerical simulations
source("plotting/plot_V_hat.R")
source("plotting/plot_PheWAS_experiment.R")
source("plotting/plot_PheWAS_experiment_BH_levels.R")
source("plotting/plot_REVIGO_experiment.R")

# plot figures for data analyses
source("plotting/plot_PheWAS_data_analysis.R")
source("plotting/plot_REVIGO_data_analysis.R")