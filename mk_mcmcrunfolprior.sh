#!/bin/bash
output_dir="/mnt/data/dayhoff/home/projects/bird_trees/bbsb/hessian/priors"

for family_directory in "$output_dir"/*/; do
    family_directory="${family_directory%/}"
    for run_number in {1..3}; do
        runFolder="$family_directory/run${run_number}"
        mkdir -p "$runFolder"

        ctlFile="$family_directory/mcmc_treearr.ctl"
        mcmcOutput="${runFolder}/run${run_number}_mcmc_pos.txt"
        output="${runFolder}/run${run_number}_out_mcmcpos.txt"

        exec >"${runFolder}/mcmctree_stdout.txt" 2>"${runFolder}/mcmctree_stderr.txt"

        sed -e "s|mcmcfile = .*|mcmcfile = ${mcmcOutput}|g" -e "s|outfile = .*|outfile = ${output}|g" "$ctlFile" > "${runFolder}/mcmc_treearr.ctl"
    done
done
