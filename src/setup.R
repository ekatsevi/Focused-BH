# source files containing functions
source("methodology/Yekutieli.R")
source("methodology/FocusedBH.R")
source("methodology/utils_DAG.R")
source("simulations/run_one_experiment.R")
source("simulations/run_one_precomputation.R")
source("utils/utils_provenance.R")
source("utils/utils_setup.R")

# load libraries
load_silently = function(packages){
  for(package in packages){
    suppressPackageStartupMessages(library(package, character.only = TRUE))    
  }
  invisible()
}
load_silently(c("readxl", "ontologyIndex", "qvalue", "cherry", 
                "structSSI", "igraph", "ape", "R.utils", "reshape2", 
                "Matrix", "multtest", "ggraph", "kableExtra",
                "janitor", "VennDiagram", "tidyverse"))

# set paths
source("utils/set_base_dir.R")
  
# download data, if not already done
system(sprintf("./download_GO_data.sh %s", base_dir))

# parse GO and UKBB data, if not already done
source("parse_GO_data.R")
source("parse_UKBB_data.R")