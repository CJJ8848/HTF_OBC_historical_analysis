#!/usr/bin/env bash
set -euo pipefail

# === INPUT: set your dump directory ===
dump_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"

# === OUTPUT files ===
out_combined="${dump_dir}/wg_kmer_mean.tsv"         # sample \t cutoff
out_modern="${dump_dir}/threshold_m57.tsv"          # starts with 'p'
out_hist="${dump_dir}/threshold_h40.tsv"            # not starting with 'p'

# Header
printf "sample\tcutoff\n" > "$out_combined"
printf "sample\tcutoff\n" > "$out_modern"
printf "sample\tcutoff\n" > "$out_hist"

shopt -s nullglob

# *_dump.txt (from jellyfish dump -c)
for f in "${dump_dir}"/*_dump.txt; do
  # skip if the glob didn't match this pattern
  [[ -e "$f" ]] || continue

  base=$(basename "$f")
  sample="${base%.*}"          # strip .txt or .tsv
  sample="${sample%_dump}"     # strip trailing _dump if present

  # Compute the filtered mean: counts > 1
  # File format: 2 columns -> kmer \t count
  mean=$(awk 'NF>=2 && $2>1 {sum+=$2; n++} END{
               if(n==0){print "NA"} else {
                 m=sum/n; printf "%.0f", (m-int(m)<0.5? int(m): int(m)+1)
               }}' "$f")

  # Write to combined
  printf "%s\t%s\n" "$sample" "$mean" >> "$out_combined"

  # Split by cohort (exactly as you asked: ONLY names starting with 'p' are modern)
  if [[ "$sample" =~ ^p ]]; then
    printf "%s\t%s\n" "$sample" "$mean" >> "$out_modern"
  else
    printf "%s\t%s\n" "$sample" "$mean" >> "$out_hist"
  fi
done

echo "✅ Wrote:"
echo "  - $out_combined"
echo "  - $out_modern"
echo "  - $out_hist"
