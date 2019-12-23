rm(list = ls())
source("setup.R")
args <- commandArgs(trailingOnly = TRUE)

# parse arguments
experiment_name = args[1]
experiment_index = as.integer(args[2])
base_dir = args[3]

# run the computation
run_one_experiment(experiment_name, experiment_index, base_dir)

cat(sprintf("Computation successfully completed!\n"))