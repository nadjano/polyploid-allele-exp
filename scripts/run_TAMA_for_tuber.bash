# Set working directory
cd /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2


# Create a directory for the Synt-specific files
# Define $Synt (ensure it is set)

# Create output directory if it doesn't exist
Synt="tuber"
mkdir -p "tuber"

tuber_ids=("SRR14995031" "SRR14995032" "SRR14995033" "SRR14995034" "SRR14995933")

# Process each BAM file
for bam_id in "${tuber_ids[@]}"; do
    bam="${bam_id}.bam"
    # Ensure file exists
    if [[ -f "$bam" ]]; then
        # Extract reads matching $Synt and save as BAM
        samtools view -h "$bam" -F 0x904 -m 500 |samtools view -b > "${Synt}/${bam%.bam}_${Synt}.bam"
    else
        echo "File $bam not found. Skipping."
    fi
done

# Define arrays for tuber and tuber IDs

#"SRR14995031" "SRR14995032" "SRR14995033" "SRR14995034" "SRR14995933")

# Define a function to process tuber and tuber BAM files
process_samples() {
  local sample_type=$1        # "tuber" or "tuber"
  local sample_ids=("${!2}")  # Array of sample IDs (tuber_ids or tuber_ids)
  
  # Merge BAM files for the given sample type
  echo "Merging BAM files for ${sample_type}..."
  echo $(printf "${Synt}/%s_${Synt}.bam " "${sample_ids[@]}")
  samtools merge -@ 8 "${Synt}/${sample_type}_${Synt}_merged.bam" $(printf "${Synt}/%s_${Synt}.bam " "${sample_ids[@]}")
  
  # Sort the merged BAM file
  echo "Sorting BAM file for ${sample_type}..."
  samtools sort -@ 8 "${Synt}/${sample_type}_${Synt}_merged.bam" -o "${Synt}/${sample_type}_${Synt}_merged_sorted.bam"
  
  # Extract only primary alignments
  echo "Extracting primary alignments for ${sample_type}..."
  samtools view -b -h -F 0x904 -m 500 "${Synt}/${sample_type}_${Synt}_merged_sorted.bam" > "${Synt}/${sample_type}_${Synt}_merged_primary.bam"
  
  # Index the primary BAM file
  echo "Indexing BAM file for ${sample_type}..."
  samtools index "${Synt}/${sample_type}_${Synt}_merged_primary.bam"
}

# Process tuber samples
process_samples "tuber" tuber_ids[@]


# We want to compare the isofroms between allels 
conda activate GSTAMA_COLLAPSE


TEST_BAM=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/${Synt}/${sample_type}_${Synt}_merged_primary.bam
# 1. Get the isofroms for each allele
FASTA=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/0_GFFREAD/Atlantic_liftoff.fasta
OUT=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_tuber

mkdir -p $OUT

TAMA_PATH=/scratch/nadjafn/git_clones/tama

#   -s S        Sorted sam/bam file (required)(if using BAM file please use -b BAM flag as well)
#   -f F        Fasta file (required)
#   -p P        Output prefix (required)
#   -x X        Capped flag: capped or no_cap

# extract only Synt_2798 reads and rerun TAMA
#samtools view -h $TEST_BAM | grep "Synt_2798_" | samtools view -b > $OUT/tuber_chr01_Synt_2798.bam

# important to remove unmapped reads
conda activate GSTAMA_COLLAPSE
python /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/scripts/tama_collapse.py -s $TEST_BAM -f $FASTA -p $OUT/tuber -b BAM -d merge_dup -i 85 -c 80 -x no_cap  -a 50 -z 50 -m 50 -e common_ends -rm low_mem

# -e longest_ends

#python $TAMA_PATH/tama_collapse.py -s $TEST_SAM  -f $FASTA -p $OUT/Synt_8877_ -x no_cap -i 60 -c 70  -vc 20 -d merge_dup -e longest_ends -m 60 -a 300 -z 150 

bedtools bed12tobed6 -i  /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test_tuber/tuber.bed  >  /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test_tuber/tuber_16.bed


# bed to gff
# conda create -n bed2gff -c bioconda bed2gff
conda activate bed2gtf
# bed2gtf --bed /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test_tuber/tuber_chr01_Synt_2798_not_long.bed --output /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test_tuber/tuber_chr01_Synt_2798_not_long.gtf --no-gene



