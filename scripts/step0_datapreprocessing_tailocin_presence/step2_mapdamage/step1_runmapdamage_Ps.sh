#!/bin/bash -l

#$ -l tmem=15G
#$ -l h_vmem=15G
#$ -l h_rt=10:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N h46_Psmapdamage 
#$ -t 1-46
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -wd /SAN/ugi/plant_genom/jiajucui/

echo "Task id is $SGE_TASK_ID"
#conda create -n  mapDamageR4
#conda activate mapdamageR4 and conda install -c conda-forge r-base=4.1 (4.2 or more is not supported) ref:https://github.com/rstudio/gt/issues/1004
#install mapDamage
#or not install the dep 5 R libs ref:https://ginolhac.github.io/mapDamage/#a1
#then use conda to install https://anaconda.org/bioconda/mapdamage2 
#conda install -c bioconda mapdamage2

source /home/jiajucui/miniconda3/bin/activate phylogeny_snp
#variables
i=$SGE_TASK_ID
samplename=$(cat samples.txt | sed -n $i'p')

wd=/SAN/ugi/plant_genom/jiajucui/phylogeny/phylogeny_snp/mapdamage/
refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta
refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
cd ${wd}/mapdamage2024_29/tops/
#echo "samplename: ${samplename}"  | tee -a ${wd}/md.log
#mapDamage -i /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2023_37HB/${samplename}*.bam -r ${refAt}
mapDamage  -i /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${samplename}.mapped_to_Pseudomonas.dd.q20.bam -r ${refPs}

