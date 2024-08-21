INPUT_DIR=/blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/out_ms/mapping_stats
SAMPLES=("Orangutan_AG06213_PAB" "RIL_updated_dm12272_01h" "Atlantic_withS")  # Add more sample names as needed

for SAMPLE in "${SAMPLES[@]}"; do
    # Concatenate files by row
    paste $INPUT_DIR/${SAMPLE}_*_P_summary.tsv > $INPUT_DIR/${SAMPLE}_P_summary.tsv
    paste $INPUT_DIR/${SAMPLE}_*_P_haplotype_summary.tsv > $INPUT_DIR/${SAMPLE}_P_haplotype_summary.tsv

    # Concatenate files by column
    cat $INPUT_DIR/${SAMPLE}_P_summary.tsv $INPUT_DIR/${SAMPLE}_P_haplotype_summary.tsv > $INPUT_DIR/${SAMPLE}_P_merged_summary.tsv
done
