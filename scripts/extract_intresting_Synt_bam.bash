# Set working directory
cd /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2

# Set the Synt identifier
#Synt="Synt_8548_"
#Synt="Synt_261_"
#Synt="Synt_5330_"
#Synt="Synt_14201_"
#Synt="Synt_30140_"
Synt="Synt_240_"
#Synt="Synt_8877_"
#Synt="Synt_18928"
# Synt="Synt_5153_"
# Synt="Synt_5408_"
# Synt="Synt_5429_"

# Create a directory for the Synt-specific files
# Define $Synt (ensure it is set)

# Create output directory if it doesn't exist
mkdir -p "$Synt"

# Process each BAM file
for bam in SRR*.bam; do
    # Ensure file exists
    if [[ -f "$bam" ]]; then
        # Extract reads matching $Synt and save as BAM
        samtools view -h "$bam" | grep "$Synt" | samtools view -b > "${Synt}/${bam%.bam}_${Synt}.bam"
        
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
leaf_ids=("SRR14993892" "SRR14993893" "SRR14993894" "SRR14993895" "SRR14996168")
tuber_ids=("SRR14995031" "SRR14995032" "SRR14995033" "SRR14995034" "SRR14995933")

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
process_samples "leaf" leaf_ids[@]

# Process tuber samples
process_samples "tuber" tuber_ids[@]

# remove the unmerged files
rm "${Synt}/SRR*_Synt_*_.bam"
# remove the merged files 
rm "${Synt}/*merged.bam" 
rm "${Synt}/*merged_sorted.bam"

# Extract sequences from FASTA file for the specific Synt identifier
awk -v Synt="$Synt" 'BEGIN {RS=">"; ORS=""} NR>1 {if ($1 ~ Synt) print ">"$0}' \
    /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/0_GFFREAD/Atlantic_liftoff.fasta \
    > "/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${Synt}.fasta"



# extract promotors
GFF=/scratch/nadjafn/Atlantic_ASE/output/liftoff_syntelogs_v6/updated_v62Atl_liftoff_a0.9s0.9_ALL.gff

# make promotor gff
grep $Synt $GFF | grep "gene" | awk -v OFS="\t" '{print $1,$2,$3,$4-1000,$4,$6,$7,$8,$9}' > /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${Synt}_promotors.gff

# chaneg gene to exon in third column
sed -i 's/gene/exon/g' /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${Synt}_promotors.gff

conda activate gffread 
# extract promotor sequences
gffread -g /scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${Synt}_promotors.gff -w  /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${Synt}_promotors.fasta
# add promotor to the fasta file names
sed -i 's/>/>promotor_/g' /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${Synt}_promotors.fasta
