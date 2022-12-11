#!/bin/bash

# genome from: /gpfs/data/sequence/cellranger-refdata/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa
# $1 path to outs folder: /gpfs/data/sequence/results/aifantislab/2020-runs/2020-03-18/cellranger/count-Audrey3/outs
# $2 out_path: /gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/data/soupercell/
# $3 sample_name: 2020-03-18-count-3
# $4 number of clusters: 5

mkdir /gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/data/soupercell/

# create file with all information needed:

outs=( $( awk -F"," 'FNR>1{print $2}' /gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/data/metadata/orig.ident_path_cluster_souporcell.csv ) )
# path to cellranger outs directory
sample=( $( awk -F"," 'FNR>1{print $1}' /gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/data/metadata/orig.ident_path_cluster_souporcell.csv ) )
# orig.ident of specific library
clusters=( $( awk -F"," 'FNR>1{print $3}' /gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/data/metadata/orig.ident_path_cluster_souporcell.csv ) )
# number of clusters to run (usually # of patients within the library)
out_path=/gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/data/soupercell/

# array length
echo ${#outs[@]}
echo ${#sample[@]}
echo ${#clusters[@]}

#run from command line

for (( i=1; i<${#outs[@]}; i++ )); do
sbatch --job-name=soupercell_"${sample[$i]}" \
    /gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/scripts/Soupercell.sh \
    ${outs[$i]} \
    $out_path \
    ${sample[$i]} \
    ${clusters[$i]}
done

