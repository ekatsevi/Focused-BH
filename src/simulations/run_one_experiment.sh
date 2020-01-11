#!/bin/bash
# Run one experiment

# load R if necessary
if [ "$(command -v module)" ]; then
  module load R
fi

script_dir=`dirname "$0"`
src_dir=$script_dir/..

cd $src_dir

# call the R script 
Rscript --vanilla call_run_one_experiment.R $1 $2 $3 $4