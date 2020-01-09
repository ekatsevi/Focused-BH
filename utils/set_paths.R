#################################################################
# SET PATHS TO BASE AND FIGURE DIRECTORIES, DEPENDING ON MACHINE
#################################################################

# machine variable can be input as a command line argument or 
# can already be present in environment 
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0){
  machine = args[1] 
} else{
  stopifnot(exists("machine"))
}

# set base directory
if(machine == "local"){
  base_dir = "/home/ekatsevi/project-files/focused-bh"  
}
if(machine == "PSC"){
  base_dir="/pylon5/ms5piap/ekatsevi/focused-bh"
}

# print base directory to standard out for bash script
if(length(args) > 0){
  cat(base_dir)  
}

# set figures directory
figures_dir = "/home/ekatsevi/Dropbox/Research/Projects/HierTest/manuscript/Reresubmission/figures"

