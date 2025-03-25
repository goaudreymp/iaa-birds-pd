#!/bin/bash 
#SBATCH --job-name=beast4
#SBATCH --output=/mnt/data/dayhoff/home/u6562250/beast/beast4.out
#SBATCH --error=/mnt/data/dayhoff/home/u6562250/beast/beast4.err
#SBATCH --partition=Standard
#SBATCH --time=240:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=u6562250@anu.edu.au
#SBATCH --mail-type=ALL

# path to your home directory 
HOME_DIR=/mnt/data/dayhoff/home/u6562250

# activate conda environment
source /opt/conda/bin/activate $HOME_DIR/.conda/envs/iaa-tree

beast  -threads 4 -prefix $HOME_DIR/beast/prioro2_ $HOME_DIR/beast/iaa_bird_09_priorsonly.xml

wait

