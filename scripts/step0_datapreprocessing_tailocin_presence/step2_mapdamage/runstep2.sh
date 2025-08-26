
for file in $(ls /path/to/mapdamageresults | grep '.mapDamage'); do
  echo $file
  Rscript step2_generateCtoTfreq_frommapdamagefolder_At.R $file $file
done


for file in $(ls /path/to/mapdamageresults/ | grep '.mapDamage'); do
  echo $file
  Rscript step2_generateCtoTfreq_frommapdamagefolder_At.R $file $file
done

