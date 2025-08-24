#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N s2025_h40
#$ -t 1-40
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

echo "Task id is $SGE_TASK_ID"
#variables
i=$SGE_TASK_ID
source ~/miniconda3/bin/activate phylogeny_snp
#variables
#need to modify: reference refPS,
#sample name and the path to store the answer
#thecollapsed file sotrage
#the path of intermediate files: ${pathtmp}
pathtmp=2025_h40_rmAT_maptoPs
answerspath=answer2025_h40_rmAT_maptoPs
mkdir -p /SAN/ugi/plant_genom/jiajucui/${answerspath}/
samplename=$(cat /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/40_OTU5.txt | sed -n $i'p')
echo "${samplename}:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
#refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta
collapsed=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.fastq.gz
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/
bwa aln -t 8 -l 1024 -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed.sai ${refAt} ${collapsed}
#Convert  reads into a standard alignment format (SAM)
bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam ${refAt} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed.sai ${collapsed}
samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_Atsam_flagstats.log 
echo "q1 What is the percentage of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
echo "reads in total (allrawreads):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
grep 'in total' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_Atsam_flagstats.log  >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


#convert mapped reads into bam with map quality score greater than 80 
samtools view -@ 2 -F 4  -q 20 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam
#rm sai, sam will be removed later
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed.sai
#for q1 What is the percentage of A.thaliana DNA?

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_At_q20_flagstats.log 
echo "reads mapped to At (-q20):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
grep 'mapped (' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_At_q20_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

#for q1.1
echo "q1.1 What is the covered genome proportion of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam &&


samtools coverage /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.bysamtoolscoverage.txt


echo "total length of At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
less /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.bysamtoolscoverage.txt | awk '{sum+=$3;} END{print sum;}' >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

echo "base mapped to At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
less /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.bysamtoolscoverage.txt | awk '{sum+=$5;} END{print sum;}' >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

#q1.2 question about the depth of At
#lenAt=119667750
echo "q1.2 What is the read depth of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
echo "read depth to At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
less /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.bysamtoolscoverage.txt | awk '{sum+=($3*$7);} END{print sum/119667750;}' >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

#then rm unsortedbam 
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam

#q2 What is the percentage of Pseudomonas DNA? after removing mapped reads to A.thaliana

####################
##separate pseudomonas from remains DNA (removal)
samtools view -bf 4 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam &&

#then rm Atsam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam


#then map to Ps

refPswithtailocin=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta
#/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/
source ~/miniconda3/bin/activate phylogeny_snp
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}
##use bwa to realign the unmapped bam file
bwa aln -t 2 -l 1024 ${refPswithtailocin} -b /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai
#_after_removal_mappedAt.bam
bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam ${refPswithtailocin} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2024_413/${samplename}_after_removal_mappedAt.bam
#Keep the mapped reads and create a compressed BAM file using samtools
samtools view -@ 2 -F 4 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam

#Sort the BAM file by chromosome and position using samtools
samtools sort -@ 2 -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.bam

samtools view -q 20 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam


#rm sai sam bam and unsort unq20 bams:
#keep for other refs
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam

samtools markdup /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_Ps_q20_flagstats.log

