#!/bin/bash
#
# Author: Gene Katsevich
# Date: 12/19/2019
# Submit jobs for Focused BH

# simulation parameters
machine="local"                 # which machine it's running on (local, ubergenno, PSC)
experiment_name="PheWAS_Fisher"    # will define input file
mode="interactive"            # interactive or batch

# set base directory depending on machine
if [ $machine == "local" ]
then
base_dir="/home/ekatsevi/project-files/focused-bh"
elif [ $machine == "PSC" ]
then
base_dir=$SCRATCH"/focused-bh"
else
  echo "Invalid machine specified!"
exit 1
fi

# paths to relevant files/directories
input_filename="input_files/input_file_"$experiment_name".R"
logs_dir=$base_dir"/logs/precomp_"$experiment_name
precomp_dir=$base_dir"/precomp/"$experiment_name
if [ ! -d "$logs_dir" ] 
then
mkdir $logs_dir
fi
if [ ! -d "$precomp_dir" ] 
then
mkdir $precomp_dir
fi

num_experiments=$(Rscript $input_filename num_precomputations)
for (( experiment_index=1; experiment_index<=$num_experiments; experiment_index++ ))
  do
echo "Submitting job for precomputation number "$experiment_index
command="./run_one_precomputation.sh $experiment_name $experiment_index $base_dir $machine"
logs_filename=$logs_dir"/"$experiment_name"_"$experiment_index".Rout"
# construct final call based on machine
if [ $machine == "local" ] 
then
if [ $mode == "batch" ] 
then
$command > $logs_filename 2> $logs_filename &
fi
if [ $mode == "interactive" ] 
then
$command
fi
fi
if [ $machine == "PSC" ]
then
if [ $mode == "batch" ] 
then
sbatch --time=02:00:00 -p RM-shared -J $experiment_index"_"$experiment_name -o $logs_filename $command
fi
if [ $mode == "interactive" ] 
then
$command
fi
fi
done