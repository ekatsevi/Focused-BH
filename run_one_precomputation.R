rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

# parse arguments
experiment_name = args[1]
b = as.integer(args[2])
base_dir = args[3]
machine = args[4]

print(machine)

source("setup.R", local = TRUE)

# run the computation
run_one_precomputation(experiment_name, b, base_dir)

cat(sprintf("Precomputation successfully completed!\n"))