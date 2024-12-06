# We want to compare the isofroms between allels 
conda activate GSTAMA_COLLAPSE

TEST_SAM=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/Synt_8877_/leaf_Synt_8877__merged_primary_mapped.sam
TEST_BAM=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/Synt_240_/tuber_Synt_240__merged_primary.bam
# 1. Get the isofroms for each allele
FASTA=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/results/01_MINIMAP2/Synt_240_/Synt_240_.fasta 
OUT=/scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test

TAMA_PATH=/scratch/nadjafn/git_clones/tama

#   -s S        Sorted sam/bam file (required)(if using BAM file please use -b BAM flag as well)
#   -f F        Fasta file (required)
#   -p P        Output prefix (required)
#   -x X        Capped flag: capped or no_cap


# important to remove unmapped reads

python /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/scripts/tama_collapse.py -s $TEST_BAM  -f $FASTA -p $OUT/Synt_240_default_no_cap_long_end -b BAM -d merge_dup -i 80 -c 75 -x no_cap -e longest_ends -a 50 -m 50 -z 50

#python $TAMA_PATH/tama_collapse.py -s $TEST_SAM  -f $FASTA -p $OUT/Synt_8877_ -x no_cap -i 60 -c 70  -vc 20 -d merge_dup -e longest_ends -m 60 -a 300 -z 150 

bedtools bed12tobed6 -i  Synt_240_.bed  >  Synt_240_6.bed   


# bed to gff
# conda create -n bed2gff -c bioconda bed2gff
conda activate bed2gtf
bed2gtf --bed /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test/Synt_240_.bed --output /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test/Synt_240_.gtf --no-gene

bed2gtf --bed /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test/Synt_240_default_no_cap_long_end.bed --output /scratch/nadjafn/Atlantic_ASE/nf-potato-ase/tama_test/Synt_240_default_no_cap_long_end.gtf --no-gene
