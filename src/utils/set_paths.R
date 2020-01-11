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
  base_dir = ".."
}

# print base directory to standard out for bash script
if(length(args) > 0){
  cat(base_dir)  
}

# set figures directory
figures_dir = "../figures"

# create directories if they don't exist
directories = c(figures_dir,
                base_dir, 
                sprintf("%s/data", base_dir), 
                sprintf("%s/data/raw", base_dir),
                sprintf("%s/data/raw/GO", base_dir), 
                sprintf("%s/data/raw/biobank", base_dir),
                sprintf("%s/data/processed", base_dir),
                sprintf("%s/data/processed/GO", base_dir),
                sprintf("%s/data/processed/biobank", base_dir),
                sprintf("%s/results", base_dir),
                sprintf("%s/precomp", base_dir),
                sprintf("%s/logs", base_dir))

for(directory in directories){
  if(!dir.exists(directory)){
    dir.create(directory)
  }
}