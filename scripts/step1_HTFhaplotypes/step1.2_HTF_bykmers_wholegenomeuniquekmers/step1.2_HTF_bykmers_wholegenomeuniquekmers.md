step0 finish, rm the shs in data done
start the ref from all to hamming, 
then run isolate jf and then hamming query isolate then run wg dist... and so on 
and continue

step0_wholegenome_7refs 
	step1_generatepaf.sh                  step2_summarypaf_corrected_byhand.sh  step3_searchanddelete.sh
#in step0 i run minimap2 to find the coordinates for the HTF region in each reference genome, then delete them to have 7 whole genome background that exclude the self HTF region for the kmer calling.

step1.1_refkmers_all_to_wgunique.sh  
#in step1.1 i run the kmer calling to output 7 sets of kmers that are unique from each other and also exclusive from the 7 representative OTU5 whole genome background.
 
step1.2_refkmers_hamming_filtering
#in step1.2 i run hamming matrix and visualized them, then run R to iteratively filter out any pairs of kmers that have distance <=1. Here we have the final kmers_unique_hamming2

step2.1_generateisolate_jf_kmersm57.sh
step2.2_generateisolate_jf_kmersh40.sh
#here i call the isolate kmers jf files
#h40 i rmdup and use the rmdup fastq to call jf

step3_wg_background_distribution
readme  step1_dump_wgdistribution_Rvisual.R  step2_calculate_mean_wg_kmercount_asthreshold.sh

#then step3 run wg distribution, first rvisual then calculate mean wg kmer count as the threshold

step4_query_and_filterby_wgthreshold
#then step4 query the kmer and filter with wg avg, run mix, observe coinfection, select the best matches...
#here i run first the summary txt including the kmer count of each HTF-isolate pair, then i calculate in R the proportion of HTF kmer to the total matched kmers to any HTF ref of one isolate, then use this as an index to rank the best matches... 
#also to plot the mix histogram of HTF and wg kmer dist, i run step1.sh to get the wg_p25.c2_depthsummary which including the kmer in long format. then run mix.R to output the dist...

step5 run coinfection, run breaking isolates and remove unsure 64, but authenticate PL0240 to be a break, then calculate HTF length group freq difference, box plot and raw data plot with time, and with geography

after cd to combine folder:
step1 combine, done
#in combine with Oantigen and now i have the final main fig
then cd ../
step2 final combine local assembly and kmer, give a confidence matrix of HTF calling

done
