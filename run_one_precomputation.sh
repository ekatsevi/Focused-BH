#!/bin/bash
# Run one experiment

if [ "$(command -v module)" ]; then
module load R
fi

Rscript --vanilla run_one_precomputation.R $1 $2 $3 $4