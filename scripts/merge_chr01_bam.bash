# Set working directory
cd /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2


# Create a directory for the Synt-specific files
# Define $Synt (ensure it is set)

# Create output directory if it doesn't exist
Synt="chr01"
mkdir -p "chr01"


# Process each BAM file
for bam in SRR*.bam; do
    # Ensure file exists
    if [[ -f "$bam" ]]; then
        # Extract reads matching $Synt and save as BAM
        samtools view -h "$bam" | grep "chr01" | samtools view -b > "${Synt}/${bam%.bam}_${Synt}.bam"
        
        # Add the sample name as RG tag
        samtools addreplacerg "${Synt}/${bam%.bam}_${Synt}.bam" \
            -r "ID:${bam%.bam}\tSM:${bam%.bam}" \
            -o "${Synt}/${bam%.bam}_${Synt}_temp.bam"
        
        # Replace the file with the updated one
        mv -f "${Synt}/${bam%.bam}_${Synt}_temp.bam" "${Synt}/${bam%.bam}_${Synt}.bam"
    else
        echo "File $bam not found. Skipping."
    fi
done

# Define arrays for leaf and tuber IDs
ids=("SRR14993892" "SRR14993893" "SRR14993894" "SRR14993895" "SRR14996168")
#"SRR14995031" "SRR14995032" "SRR14995033" "SRR14995034" "SRR14995933")

# Define a function to process leaf and tuber BAM files
process_samples() {
  local sample_type=$1        # "leaf" or "tuber"
  local sample_ids=("${!2}")  # Array of sample IDs (leaf_ids or tuber_ids)
  
  # Merge BAM files for the given sample type
  echo "Merging BAM files for ${sample_type}..."
  samtools merge -@ 8 "${Synt}/${sample_type}_${Synt}_merged.bam" $(printf "${Synt}/%s_${Synt}.bam " "${sample_ids[@]}")
  
  # Sort the merged BAM file
  echo "Sorting BAM file for ${sample_type}..."
  samtools sort -@ 8 "${Synt}/${sample_type}_${Synt}_merged.bam" -o "${Synt}/${sample_type}_${Synt}_merged_sorted.bam"
  
  # Extract only primary alignments
  echo "Extracting primary alignments for ${sample_type}..."
  samtools view -b -h -F 2308 -m 300 "${Synt}/${sample_type}_${Synt}_merged_sorted.bam" > "${Synt}/${sample_type}_${Synt}_merged_primary.bam"
  
  # Index the primary BAM file
  echo "Indexing BAM file for ${sample_type}..."
  samtools index "${Synt}/${sample_type}_${Synt}_merged_primary.bam"
}

# Process leaf samples
process_samples "_chr01" ids[@]

