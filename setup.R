# paths
if(machine == "local"){
  base_dir = "/home/ekatsevi/project-files/focused-bh"  
}
if(machine == "PSC"){
  base_dir="/pylon5/ms5piap/ekatsevi/focused-bh"
}

# load libraries
load_silently = function(packages){
  for(package in packages){
    suppressPackageStartupMessages(library(package, character.only = TRUE))    
  }
  invisible()
}

load_silently(c("readxl", "ontologyIndex", "qvalue", "cherry", 
                "structSSI", "igraph", "ape", "R.utils", "reshape2", 
                "tidyverse"))


# download data and run preprocessing
system(sprintf("./download_GO_data.sh %s", base_dir))
source("parse_GO_files.R")
# source("parse_ICD_files.R")

# auxiliary functions
source("aux_DAG.R")
source("simulation_utils.R")
source("utils_provenance.R")