#!/bin/bash

# path to your home directory
HOME_DIR=/mnt/data/dayhoff/home/u6562250

# activate conda environment
source /opt/conda/bin/activate $HOME_DIR/.conda/envs/iaa-tree

# list folder
fasta="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/alignments"
tree="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/trees"

output="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/baseml"
template_ctl="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/mcmc_baseml.ctl"

for fasta_file in $fasta/*.fasta;
do
family=$(basename $fasta_file .fasta)

echo "$family"

family_directory=$output/$family
mkdir -p $family_directory

cp $template_ctl $family_directory

sed -i "s|seqfile = fullpass_concat.fasta|seqfile = $fasta/$family.fasta|" "$family_directory/mcmc_baseml.ctl"
sed -i "s|treefile = fullcalib_oliveros.tre|treefile = $tree/$family.newick|" "$family_directory/mcmc_baseml.ctl"

    echo "Created directory for $family_name"

done


