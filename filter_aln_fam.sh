#!/bin/bash

# path to your home directory
HOME_DIR=/mnt/data/dayhoff/home/u6562250

# activate conda environment
source /opt/conda/bin/activate $HOME_DIR/.conda/envs/iaa-tree

# list folder
list="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/list"
fasta="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/alignments/concat_trimmed_090.fasta"
output="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/alignments"

	for list_file in $list/*.txt;
	do
		list_filename=$(basename $list_file)

		prefix=${list_filename%List*}

		output_filename=${prefix}.fasta
		output_path=$output/$output_filename

	    	faSomeRecords $fasta $list_file $output_path

		echo "Subsetted fasta created for $list_filename" 
	done
done

