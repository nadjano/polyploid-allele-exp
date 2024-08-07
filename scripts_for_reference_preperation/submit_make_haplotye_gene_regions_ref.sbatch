#!/bin/bash
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --job-name=syntelogs_class
#SBATCH --output=syntelogs.out
#SBATCH --error=/syntelogs.err
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=50gb

module load python/3.11
module load samtools/1.19
module load gffread/0.12.7 

### Set Directories
PROJ="directory of current project"
SCRIPTS=$PROJ/scripts/construct_haplotype_generegions_and_compare
OUTD=$PROJ/transcripts_classification
ROZ=$OUTD/roz
REF_DIR=$PROJ/spuddb_reference_data

mkdir -p $ROZ


SYNTELOGS=$REF_DIR/syntelogs_atl0to4.csv
TRANSCRIPTS_GFF=$REF_DIR/ATL_v3.hc_gene_models.repr.gff3.gz
GENOME_FASTA=$REF_DIR/ATL_v3.asm.fa # Note gunzip and then bgzip
# Filter Syntelogs file to only get atlantic gene ids
grep "Soltu.Atl_v3" $SYNTELOGS > $ROZ/syntelogs_atl0to4_only_atl.csv


# We want for every gene get the number of allelic variants and copies and SNPs between them
# We dont have unified gene ids accross the haplotyes so we use syntelogs ids as "new gene ids" (e.g Synt_1_hap_x[number of allels in group])

python $SCRIPTS/rename_geneIDs_cg.py \
                        --syntelogs_file $ROZ/syntelogs_atl0to4_only_atl.csv \
                        --gff3 $TRANSCRIPTS_GFF \
                        --gff3_out $REF_DIR/ATL_v3.hc_gene_models_syntIDs.repr.gff3


# Index fasta to speed up gffread extraction
samtools faidx $REF_DIR/ATL_v3.asm.fa

# Extract mRNA rows only
grep "mRNA" $REF_DIR/ATL_v3.hc_gene_models_syntIDs.repr.gff3 > $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr.gff3
# Rename mRNA to exon for gffread gene region extraction
sed -i 's/mRNA/exon/g' $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr.gff3

# extract gene regions with adapted gff3 file
gffread $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr.gff3 -g $GENOME_FASTA -v -E -w $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr.fa --w-add 150

# Put all sequence per gene region in one line to make is easier to extract sequences
perl -pe '/^>/ ? print "\n" : chomp' $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr.fa | tail -n +2 > $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr_one_line_allhap.fa


# Extract all sequences that have _1, or _2, _3, _4, _0 in the sequence name in a loop
for i in {0..4}
do
    grep "_${i}G" -A1 $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr_one_line_allhap.fa > $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr_one_line_hap$i.fa
done
# Add the scaffold sequences to the haplotype 1 sequences
grep "Soltu.Atl_v3.S" -A1 $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr_one_line_allhap.fa >> $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr_one_line_hap1.fa

# copy to hap_genomes dir for nf-ASE-mapping-comparison
cp $REF_DIR/ATL_v3.hc_gene_models_only_mRNA_syntIDs.repr_one_line_*.fa $REF_DIR/hap_genomes