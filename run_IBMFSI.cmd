#!/bin/bash
#PBS -N "job_case_name"
#PBS -e errorfile.err
#PBS -o logfile.log
#PBS -q gpu
#PBS -l select=1:host=node8:ncpus=2:ngpus=1
#PBS -V
cd $PBS_O_WORKDIR


ulimit -s unlimited

export CUDA_VISIBLE_DEVICES=MIG-983a1c93-e5cd-5b44-8b3c-ca756e58135c
time ./make.out > output.txt 
