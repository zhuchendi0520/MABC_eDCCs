#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 3-00:00:00
#SBATCH --mem=30G

#tree builder to iqtree

input_fa=$1
output_name=$(basename $input_fa)
output_name=${output_name%_core_gene_alignment.aln}

#run_gubbins.py --prefix $output_name --tree-builder raxml $input_fa \
#    --threads 24

run_gubbins.py --prefix $output_name --tree-builder iqtree $input_fa \
    --threads 24