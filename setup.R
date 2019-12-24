# paths
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
library(readxl, quietly = TRUE)
library(ontologyIndex, quietly = TRUE)
library(qvalue, quietly = TRUE)
library(cherry, quietly = TRUE)      # to run the Structured Holm procedure (Meijer and Goeman, 2015b)
library(structSSI, quietly = TRUE)   # to perform Yekutieli test
library(igraph, quietly = TRUE)      # to perform Yekutieli test
library(ape, quietly = TRUE)         # to perform Yekutieli test
library(R.utils, quietly = TRUE)
library(reshape2, quietly = TRUE)
suppressPackageStartupMessages(library(tidyverse))

# auxiliary functions
source("aux_DAG.R")
source("simulation_utils.R")