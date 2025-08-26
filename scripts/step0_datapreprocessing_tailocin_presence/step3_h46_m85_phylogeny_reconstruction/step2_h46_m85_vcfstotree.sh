#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -N all131olden 
#$ -wd /path/
#$ -V
#$ -e /path/logs/
#$ -o /path/logs/
#wd=/path/phylogeny/phylogeny_snp
wd=/path/phylogeny/phylogeny_snp/phy/anew2024/vcfandphy/all131/
cd ${wd}
output_file=${wd}/SNPstats.txt
samplenumber=131

source ~/miniconda3/bin/activate r2t

#merge remove the 7 temporarry
bcftools merge /path/phylogeny/phylogeny_snp/vcfs/samples46_85_aftertrimnsqs10/*filtered.sorted.vcf.gz  -Oz -o ${wd}/oldmergesnp_filtered.anew2024_all131.vcf.gz

zless ${wd}/oldmergesnp_filtered.anew2024_all131.vcf.gz | sed 's#/path/4_mapping_to_pseudomonas/trimnstrimqualities10/##g' | sed 's/_removalAt_mapped_to_ps.sorted.q20.bam//g' | gzip > ${wd}/mergesnp_filtered.anew2024_all131.vcf.gz
 
rm ${wd}/oldmergesnp_filtered.anew2024_all131.vcf.gz

bcftools index -t ${wd}/mergesnp_filtered.anew2024_all131.vcf.gz
# the sites after completeness filtering

nofall=$(bcftools view -g ^miss -H ${wd}/mergesnp_filtered.anew2024_all131.vcf.gz | wc -l)

#python /path/phylogeny/phylogeny_snp/vcf2phylip/vcf2phylip.py -i ${wd}/mergesnp_filtered.anew2024_all131.vcf.gz -m 85

#to know the nomulti site number
# bgzip nomulti for invar count:
bcftools view -M2 ${wd}/mergesnp_filtered.anew2024_all131.vcf.gz | bgzip > ${wd}/nomultisnp_filtered.anew2024_all131.vcf.bgzip.gz
bcftools index -t ${wd}/nomultisnp_filtered.anew2024_all131.vcf.bgzip.gz && echo 'bgzip'

nofnomulti=$(bcftools view -g ^miss -H ${wd}/nomultisnp_filtered.anew2024_all131.vcf.bgzip.gz | wc -l)
#biallelic

bcftools view -m2 -M2 ${wd}/mergesnp_filtered.anew2024_all131.vcf.gz | gzip > ${wd}/bialleliconly_filtered.anew2024_all131.vcf.gz

nofbiallelic=$(bcftools view -g ^miss -H ${wd}/bialleliconly_filtered.anew2024_all131.vcf.gz | wc -l)

python /path/phylogeny/phylogeny_snp/vcf2phylip/vcf2phylip.py -i ${wd}/bialleliconly_filtered.anew2024_all131.vcf.gz -m ${samplenumber}

[ ! -f $output_file ] && echo 'allmergedSNPs,nomultiSNPs,biallelicSNPs' > $output_file 
echo $nofall,$nofnomulti,$nofbiallelic >> $output_file

#the cor of each 250213 SNPs
bcftools view -g ^miss ${wd}/bialleliconly_filtered.anew2024_all131.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > ${wd}/fullinfo_bialleliccor.txt

cd ${wd}

iqtree -T AUTO -s ${wd}/bialleliconly_filtered.anew2024_all131.min131.phy  -m MFP -B 1000 -alrt 1000 | tee -a ${wd}/phy/anew2024/vcfandphy/all131/iqtree.log

