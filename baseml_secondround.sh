#!/bin/bash
#SBATCH --job-name=baseml_array
#SBATCH --output=array_%A_%a.out
#SBATCH --array=1-77
#SBATCH --cpus-per-task=2
#SBATCH --time=60-00:00:00
#SBATCH --mail-user=u6562250@anu.edu.au
#SBATCH --mail-type=ALL

desired_folders=("ALAUDIDAE" "ATRICHORNITHIDAE" "CORVIDAE" "EMBERIZIDAE" "FRINGILLIDAE" "IRENIDAE" "LANIIDAE" "MENURIDAE" "PARADISAEIDAE" "PARIDAE" "PHYLLOSCOPIDAE" "PTILONORHYNCHIDAE" "SITTIDAE")

directories=("/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/baseml"/*)

current_folder="${directories[$SLURM_ARRAY_TASK_ID - 1]}"

cd "$current_folder"

echo "Checking in $current_folder"

# Check if the current folder is in the desired list
folder_name=$(basename "$current_folder")
if [[ " ${desired_folders[@]} " =~ " $folder_name " ]]; then
    echo "Running baseml for $ctl_file..."
    ctl_file="$current_folder/tmp0001.ctl"
    if [ -e "$ctl_file" ]; then
        /mnt/data/dayhoff/home/u6562250/bin/baseml "$ctl_file"
    else
        echo "No tmp0001.ctl found in $current_folder."
    fi
else
    echo "Skipping $current_folder."
fi
