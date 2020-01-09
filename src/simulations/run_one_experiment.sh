#!/bin/bash
# Run one experiment

# load R if necessary
if [ "$(command -v module)" ]; then
  module load R
fi

# call the R script 
Rscript --vanilla call_run_one_experiment.R $1 $2 $3 $4