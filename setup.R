# paths
print(machine)

if(machine == "local"){
  base_dir = "/home/ekatsevi/project-files/focused-bh"  
}
if(machine == "PSC"){
  base_dir="/pylon5/ms5piap/ekatsevi/focused-bh"
}

# download data and run preprocessing
system(sprintf("./download_GO_data.sh %s", base_dir))
source("parse_GO_files.R")

# libraries
library(readxl)
library(ontologyIndex)
library(qvalue)
library(cherry)      # to run the Structured Holm procedure (Meijer and Goeman, 2015b)
library(structSSI)   # to perform Yekutieli test
library(igraph)      # to perform Yekutieli test
library(ape)         # to perform Yekutieli test
library(R.utils)
library(reshape2)
library(tidyverse)

# auxiliary functions
source("aux_DAG.R")
source("simulation_utils.R")