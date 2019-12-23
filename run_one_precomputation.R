rm(list = ls())
source("setup.R")
args <- commandArgs(trailingOnly = TRUE)

# parse arguments
experiment_name = args[1]
b = as.integer(args[2])
base_dir = args[3]

# run the computation
run_one_precomputation(experiment_name, b, base_dir)

cat(sprintf("Precomputation successfully completed!\n"))