#!/bin/bash 
#SBATCH --job-name=beast2
#SBATCH --output=/mnt/data/dayhoff/home/u6562250/beast/beast2.out
#SBATCH --error=/mnt/data/dayhoff/home/u6562250/beast/beast2.err
#SBATCH --partition=Standard
#SBATCH --exclude=wright,fisher
#SBATCH --time=480:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=u6562250@anu.edu.au
#SBATCH --mail-type=ALL

# path to your home directory 
HOME_DIR=/mnt/data/dayhoff/home/u6562250

# activate conda environment
source /opt/conda/bin/activate $HOME_DIR/.conda/envs/iaa-tree

beast  -resume -threads 4 -prefix prior2_ $HOME_DIR/beast/iaa_bird_09_priors.xml

wait

