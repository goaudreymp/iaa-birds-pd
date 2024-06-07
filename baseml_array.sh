#!/bin/bash
#SBATCH --job-name=baseml_array
#SBATCH --output=array_%A_%a.out
#SBATCH --array=1-77
#SBATCH --cpus-per-task=2
#SBATCH --time=60-00:00:00
#SBATCH --mail-user=u6562250@anu.edu.au
#SBATCH --mail-type=ALL

directories=("/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/baseml"/*)

current_folder="${directories[$SLURM_ARRAY_TASK_ID - 1]}"

cd $current_folder

echo "Checking in $current_folder"

# Check if tmp001.ctl exists in the current folder
ctl_file="$current_folder/tmp0001.ctl"
if [ -e "$ctl_file" ]; then
    echo "Running baseml for $ctl_file..."
    /mnt/data/dayhoff/home/u6562250/bin/baseml "$ctl_file"
else
    echo "No tmp0001.ctl found in $current_folder."
fi
