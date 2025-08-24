by local assembly:

historical: later follow before the shtailocin extract read ids, then assembly then map using minimap2
old wd: /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/step1 and step2
now we dont need the seq results, stop at paf, then using the same pipe like in modern57 to summary...

now we have 9 recovered...


Before, for the HTF extraction, we simply assembled the historical tailocin regions using raw reads mapped to the reference tailocin with seven HTF haplotypes and then minimap the contigs to reference to see if there is any HTF haplotypes present in particular historical sample. 

First step we successfully assembled 30 out of 40 historical samples. 

Second step, we have 24 historical HTFs present in the assemblies. 

Third step, we have 10 good quality historical HTFs survived a 65% covered proportion threshold to their certain reference haplotypes. ((and 57/74 isolates have a covered proportion to the haplotype reference of certain HTF higher than 65%...) so all 47 modern samples are good)

And the last step, i translated the 10 nucletide seqs to amino acid and found one PL0046 is incomplete, thus deleted it. Finally we have 9 historical HTFs, with high confidence.

Later, check Kmer accuracy.. Kmer we have 36 recovered...





modern: run 57 with minimap2 assembly against query... then we have the HTF haplotype by local assembly, select the top matched HTF and TFA by covered_prop (we can see most of them have >98% covered prop to the top HTF ref, along with a good presence of TFA of the same ref, even though sometimes the best TFA is not the best HTF ref).
later check confidence matrix with kmer.
