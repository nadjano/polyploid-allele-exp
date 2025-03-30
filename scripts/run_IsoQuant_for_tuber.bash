# Set working directory
FASTQ=/scratch/nadjafn/Input_data/atlantic_ONT
OUT=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/isoquant_tuber
FASTA=/scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa
GTF=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/isoquant_tuber/updated_v62Atl_liftoff_a0.9s0.9_ALL.corrected.gtf
sample_ids=("SRR14995031" "SRR14995032" "SRR14995033" "SRR14995034" "SRR14995933")

conda activate isoquant
isoquant.py --reference $FASTA \
  --fastq $(printf "${FASTQ}/%s.fastq " "${sample_ids[@]}") \
  --data_type nanopore -o $OUT --genedb $GTF --labels tuber tuber tuber tuber tuber



