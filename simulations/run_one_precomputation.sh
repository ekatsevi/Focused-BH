#!/bin/bash
# Run one precomputation

if [ "$(command -v module)" ]; then
module load R
fi

Rscript --vanilla call_run_one_precomputation.R $1 $2 $3 $4