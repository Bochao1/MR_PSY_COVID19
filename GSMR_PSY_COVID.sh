#!/bin/bash
#SBATCH -t 05:30:00
#SBATCH --mem=300G
#SBATCH -o log.out
#SBATCH -e errlog.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=b.lin@umcutrecht.nl

dir=/hpc/hers_en/blin/Publicdata/COVID
gcta=/hpc/hers_en/blin/program/gcta_1.93.0beta/gcta64

###reference files, 5000 random selected UKB samples plink file list 
plinklist=/hpc/hers_en/blin/Publicdata/ENIGMA/GSMR/reference/ref.file.list
PSY_list=$dir/pgc/psy.list
COVID_list=$dir/exposure/covid.list

##perfom bidrectional forward univariable MR

$gcta \
--mbfile  $(plinklist) \
--gsmr-file ${PSY_list} ${COVID_list} \
--gsmr-snp-min 2 \
--clump-r2 0.01 \
--effect-plot \
--gsmr-direction 0 \
--out $dir/Uni_MR_PGC_COVID \
--diff-freq 0.2 \
--thread-num 10
 