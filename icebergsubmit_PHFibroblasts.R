#!/bin/bash
#$ -pe openmp 8
#memory requests are per-core
#$ -l rmem=8G -l mem=8G
#Prefer the hidelab queue but spill over to over queues if it is full
#$ -P hidelab
#$ -j y

module load apps/gcc/5.2/bcbio/0.9.6a
work_dir='/shared/hidelab2/user/mdp15cmg/TDP-43/PH_Fibroblasts'

#Seq.Reads file directories
r1_files=$work_dir/input/Read1
r2_files=$work_dir/input/Read2

#Read in seq reads
r1=($(find $r1_files -type f -name "*.gz"|sort -n))
r2=($(find $r2_files -type f -name "*.gz"|sort -n))

#Download the best-practice template file for RNAseq experiment
echo "DOWNLOADING TEMPLATE"
bcbio_nextgen.py -w template illumina-rnaseq PH_bcbio

#Edit the template
echo "EDITTING TEMPLATE"
#Switch to using star
sed -i 's/tophat2/star/g' $work_dir/PH_bcbio/config/PH_bcbio-template.yaml

#Initialise the main analysis
echo "INITIALISING ANALYSIS"
bcbio_nextgen.py -w template $work_dir/PH_bcbio/config/PH_bcbio-template.yaml $work_dir/PH_bcbio/config/PH_bcbio.csv ${r1[@]} ${r2[@]}

#Perform the analysis
echo "PERFOMING ANALYSIS"
cd $work_dir/PH_bcbio/work
bcbio_nextgen.py -n 8 ../config/PH_bcbio.yaml

