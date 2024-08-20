#!/bin/bash
# list folder
fasta="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/alignments/dummy"

for fasta_file in $fasta/*.fasta;
do
family=$(basename $fasta_file .fasta)

tree="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/trees/calib"
output="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/priors"
template_ctl="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/mcmc_treearr.ctl"
baseml="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/baseml"

family_directory=$output/$family
mkdir -p $family_directory

cp $template_ctl $family_directory

sed -i "s|seqfile = fullpass_concat.fasta|seqfile = $fasta/$family.fasta|" "$family_directory/mcmc_treearr.ctl"
sed -i "s|treefile = fullcalib_oliveros.tre|treefile = $tree/$family.newick|" "$family_directory/mcmc_treearr.ctl"
sed -i "s|usedata = 2 /in.BV|usedata = 0|" "$family_directory/mcmc_treearr.ctl"

    echo "Created directory for $family"

done
