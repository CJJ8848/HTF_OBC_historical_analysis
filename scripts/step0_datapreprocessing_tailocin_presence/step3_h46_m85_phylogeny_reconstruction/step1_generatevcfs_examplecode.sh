#!/bin/bash -l

#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=20:00:0
#$ -wd /path/
#$ -V
#$ -e /path/logs/
#$ -o /path/logs/

#variables
i=${1}
samplename=$(cat /path/phylogeny/samplelist/all131.txt | sed -n $i'p')
wd=/path/phylogeny/phylogeny_snp
echo "${samplename}:" | tee -a ${wd}/snp.log
refPs=/path/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta

#from bam to snp calling
#-Ou output bcf files
#-c consensus model (not -m multiallels model since we have a dominant single strain with few allel)
#-Oz is the compressed gz file
#we need a MSA so donot use -v only show variant sites only
#view -i '%QUAL>=20%' phred score filtration, for 1 sample it is exactly same a s -q20, just double check
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools mpileup -Ou -f ${refPs} /path/4_mapping_to_pseudomonas/all1315/${samplename}.mapped_to_Pseudomonas.dd.q20.bam | /SAN/ugi/plant_genom/software/bcftools-1.11/bcftools call -c --skip-variants indels -Oz --ploidy 1 | /SAN/ugi/plant_genom/software/bcftools-1.11/bcftools view -i '%QUAL>=20' -Oz > ${wd}/vcfs/groups/all131/${samplename}.ploidy1_filtered.vcf.gz | tee -a ${wd}/snponly.log &&
echo 'vcf generated'

/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools sort /${wd}/vcfs/groups/all131/${samplename}.ploidy1_filtered.vcf.gz -Oz -o ${wd}/vcfs/groups/all131/${samplename}.ploidy1_filtered.sorted.vcf.gz| tee -a ${wd}/snponly.log &&
echo 'vcf sort' 

#index with format tbi (instead of csi)
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -t ${wd}/vcfs/groups/all131/${samplename}.ploidy1_filtered.sorted.vcf.gz| tee -a ${wd}/snponly.log &&
echo 'vcf index tbi'

rm  ${wd}/vcfs/groups/all131/${samplename}.ploidy1_filtered.vcf.gz 
