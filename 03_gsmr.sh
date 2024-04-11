#!/bin/bash

#PBS -N gsmr
#PBS -q normal
#PBS -l select=2:ncpus=24:mem=96GB
#PBS -o log/
#PBS -e log/
#PBS -l walltime=2:00:00
#PBS -P ProjectID


##################
### Parameters ###
##################
working="path_to_working_directory"

bfile="${working}/1KG_Phase3_refpanel/1KG_all_phase3_ns"
exposure="${working}/gsmr_exposure_metabolite.txt"
outcome="${working}/gsmr_outcome_cad.txt"
out="${working}/metabolite_cad_gsmr"

##################
### GSMR ###
##################
gcta64 --bfile ${bfile} --gsmr-file ${exposure} ${outcome} --gsmr-direction 0 --diff-freq 0.5 --gwas-thresh 1e-05 --clump-r2 0.05 --gsmr-snp-min 10 --heidi-thresh 0.01 --out ${out} > ${out}.log 2>&1
