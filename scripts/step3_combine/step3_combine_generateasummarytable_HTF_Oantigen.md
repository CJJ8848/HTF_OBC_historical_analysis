#the rest of five genes:
binary_otherfivegenes=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical46_modern57_bybamcov/final_binary_matrix.tsv

#first generate a summary table including sample HTFgroup_Oantigen_PA wfgD_PA rmlC1_PA tagG1_PA tagH1_PA spsA_PA espE2_PA espE2length note (input manually afterwards):
# HTF group are from: first col is sample name, second col is HTF haplotypes, also the sample list is the samples from this file
/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/form57_h37_HTF_bykmer_tmp_addthe2modernbymapping.txt
#copy from /Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step3_querynewref_distribution/wgfinal_combinedepthandprop_avg/s4_epse2info/HTF_bykmer_final_filtered_corrected64.tsv

#HTFgroup_Oantigen_PA: if HTF is *(1830) and *(1383), then the value is OBC_absent, if is *(1803) or *(1245), the value is OBC_present
#wfgD_PA rmlC1_PA tagG1_PA tagH1_PA spsA_PA these can be found from binary_otherfivegenes=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical46_modern57_bybamcov/final_binary_matrix.tsv row is genes, col are samples... if no sample then type NA

#espE2_PA should check this /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt, if there is a seq, then is present, set to 1 if not 0.
#espE2length should be counted from the same file completeopenframe.18empty_39fasta.txt, keep the note, would be type manually...

output to  /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/combined_HTF_oantigen_m57_h37_needOantigenfivegenePAin2modern.txt


#then appendx one row named epsE2 and if espE2 is present (39), set a 1 to the binary matrix.
binary_otherfivegenes=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical46_modern57_bybamcov/final_binary_matrix.tsv
output as /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/final_sixgenes_binary_matrix.tsv
#ls /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/
#39fasta_emboss_transeq-I20250724-161752-0560-416209-p1m.out.txt  completeopenframe.18empty_39fasta.txt  readme  summarytable

WORKDIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy" 
