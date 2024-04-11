#!/bin/bash

#PBS -N smr
#PBS -q normal
#PBS -l select=1:ncpus=24:ompthreads=24:mem=96GB
#PBS -o log/
#PBS -e log/
#PBS -l walltime=12:00:00
#PBS -P ProjectID

##################
### Parameters ###
##################
working="path_to_working_directory"

bfile="${working}/1KG_Phase3_refpanel/1KG_all_phase3_ns"
sum_stat="${working}/metabolite_sumstat.txt"
eqtl="${working}/GTEx_WholeBlood_ciseQTL/Whole_Blood"
out="${working}/metabolite_gtex_WBciseqtl_smr"

##################
### SMR ###
##################
smr_Linux --bfile ${bfile} --maf 0.01 --diff-freq-prop 0.25 --gwas-summary ${sum_stat} --beqtl-summary ${eqtl} --smr-multi --thread-num 24 --out ${out} > ${out}.smr.log 2>&1
