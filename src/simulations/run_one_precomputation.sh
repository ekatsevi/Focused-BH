#!/bin/bash
# Run one precomputation

if [ "$(command -v module)" ]; then
module load R
fi

script_dir=`dirname "$0"`
src_dir=$script_dir/..

cd $src_dir

Rscript --vanilla call_run_one_precomputation.R $1 $2 $3 $4