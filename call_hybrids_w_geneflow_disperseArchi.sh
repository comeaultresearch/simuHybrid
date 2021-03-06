#!/bin/bash
#SBATCH --job-name=call_hybrids_w_geneflow_disperse
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH -N 1                          # Esure that all cores are on one machine
#SBATCH --time=6-12:00                # Run time
#SBATCH --mem=24g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /proj/matutelb/users/aacomeault/09_logs/call_hybrids_w_geneflow_disperse-%A.out
#SBATCH -e /proj/matutelb/users/aacomeault/09_logs/call_hybrids_w_geneflow_disperse-%A.err

module add python/2.7.12 

my_s_dmi=$1
my_s_add=$2
my_h=$3
my_r=$4
my_n=$5
my_m=$6

python /proj/matutelb/users/aacomeault/00_stripts/simupop/simuHybrid_hybrids_w_geneflow_disperseArchi_NEW.py $my_s_dmi $my_s_add $my_h $my_r $my_n $my_m \
>  /proj/matutelb/users/aacomeault/02_out/simupop/simuHyb_geneflow/disperseArch_sel_1/'simuHyb_geneflow_disperseArch_sel1_pmix_data_n'$my_n'_m'$my_m'_r'$my_r'_sDMI'$my_s_dmi'_sADD'$my_s_add'.txt'


