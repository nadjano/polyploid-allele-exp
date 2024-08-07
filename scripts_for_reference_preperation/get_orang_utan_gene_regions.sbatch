#!/bin/bash
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --job-name=syntelogs_class
#SBATCH --output=orang_utan.out
#SBATCH --error=orang_utan.err
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=50gb


RPOJ="nf-ASE-mapping-comparison directory"

module load sra
module load gffread
module load samtools

mkdir -p $RPOJ/orang_utan_references
cd $RPOJ/orang_utan_references

# remove previous fai indices 
rm *.fai

wget https://eichlerlab.gs.washington.edu/public/primate_t2t_liftoff_refseq/AG06213_PAB.hap1.refseq.liftoff.v2.gff

wget https://eichlerlab.gs.washington.edu/public/primate_t2t_liftoff_refseq/AG06213_PAB.hap2.refseq.liftoff.v2.gff

# get the ncbi gff for the mitochondrial genes 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gff.gz

zcat GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gff.gz | grep "NC_002083.1" > mitochrondrial_genes.gff


# primary reference
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/685/GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt/GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt_genomic.fna.gz

# alt reference
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz


# fasterq-dump SRR27438208

# extract all lines with "chrY_hap2_hsaY" from AG06213_PAB.hap2.refseq.liftoff.v2.gff and append it to AG06213_PAB.hap1.refseq.liftoff.v2.gff
cp AG06213_PAB.hap1.refseq.liftoff.v2.gff AG06213_PAB.hap1.refseq.liftoff.v2_withY.gff
grep "chrY_hap2_hsaY" AG06213_PAB.hap2.refseq.liftoff.v2.gff >> AG06213_PAB.hap1.refseq.liftoff.v2_withY.gff
# add the mitochondrial genes to the gff file
cat mitochrondrial_genes.gff >> AG06213_PAB.hap1.refseq.liftoff.v2_withY.gff

# Extract only the mRNA rows
grep "ID=gene" AG06213_PAB.hap1.refseq.liftoff.v2_withY.gff > AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff

grep "ID=gene" AG06213_PAB.hap2.refseq.liftoff.v2.gff > AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff

# replace " " with "_" in the description column
sed -i 's/ /_/g' AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff
sed -i 's/ /_/g' AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff

# get chromosome names
grep -o -E "^[a-zA-Z0-9_\.]+" AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff  | uniq  > hap1_gff_chromosomes.txt

grep -o -E "^[a-zA-Z0-9_\.]+" AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff  | uniq  > hap2_gff_chromosomes.txt

grep   -E "^>CM[0-9]+\.[0-9]" GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt_genomic.fna  > hap2_fasta_chromosomes.txt 

grep   -E ">NC_[0-9]+\.[0-9]" GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna > hap1_fasta_chromosomes.txt 

sed -i '/chr18_hap1_hsa16_ra/d'  hap1_gff_chromosomes.txt

paste hap2_fasta_chromosomes.txt hap2_gff_chromosomes.txt > hap2_chromosomes.txt

paste hap1_fasta_chromosomes.txt hap1_gff_chromosomes.txt > hap1_chromosomes.txt

module load python
python $PROJ/scripts_for_reference_preperation/rename_chromosomes_inFasta.py hap2_chromosomes.txt GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt_genomic.fna GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt_genomic_renamed.fna

python $PROJ/scripts_for_reference_preperation/rename_chromosomes_inFasta.py hap1_chromosomes.txt GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic_renamed.fna

# Replace gene in the thrid column with exon
awk 'BEGIN{OFS="\t"} {if ($3 == "gene") $3 = "exon"; print}' AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff > temp && mv temp AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff

awk 'BEGIN{OFS="\t"} {if ($3 == "gene") $3 = "exon"; print}'  AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff > temp && mv temp AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff

# Add hap prefix to gene names
sed -i 's/ID=gene/ID=hap1_gene/g' AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff
sed -i 's/ID=gene/ID=hap2_gene/g' AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff

# remove dupliated rows from gff files
awk '!seen[$0]++' AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff > temp && mv temp AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff
awk '!seen[$0]++' AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff > temp && mv temp AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff


cat AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff > AG06213_PAB.hap1.hap2.refseq.liftoff.v2_only_gene.gff

# remove genes in gff that are not in fasta
grep -v "chr18_hap1_hsa16_random" AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff > temp && mv temp AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff

grep -v "chrY_hap2_hsaY" AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff > temp && mv temp AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff

# Extract gene regions
gffread AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff -g GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic_renamed.fna  -E -w AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.fa --w-add 150

# Extract gene regions
gffread AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.gff -g GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt_genomic_renamed.fna -E -w AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.fa --w-add 150


# remove the chrX and chrY and mito from the gff files
grep -E -v  "chrX|chrY|"NC_002083.1"" AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.gff >  AG06213_PAB.hap1.refseq.liftoff.v2_only_gene_noXYM.gff

# Extract gene regions
gffread AG06213_PAB.hap1.refseq.liftoff.v2_only_gene_noXYM.gff -g GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic_renamed.fna  -E -w AG06213_PAB.hap1.refseq.liftoff.v2_only_gene_noXYM.fa --w-add 150

# concat the hap1 and hap2 gene regions
cat AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.fa AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.fa > AG06213_PAB.hap1.hap2.refseq.liftoff.v2_only_gene.fa

cat AG06213_PAB.hap1.refseq.liftoff.v2_only_gene_noXYM.fa AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.fa > AG06213_PAB.hap1.hap2.refseq.liftoff.v2_only_gene_noXYM.fa

cp *only_gene.fa $RPOJ/nf-ASE-mapping-comparison/hap_genomes

#cp SRR27438208.fastq $RPOJ/nf-ASE-mapping-comparison/fastq_reads


# get the shared genes
grep ">" AG06213_PAB.hap2.refseq.liftoff.v2_only_gene.fa | sort | cut -d'_' -f2 > AG06213_PAB.hap2.genes.tsv

grep ">" AG06213_PAB.hap1.refseq.liftoff.v2_only_gene.fa | sort | cut -d'_' -f2 > AG06213_PAB.hap1.genes.tsv

sort AG06213_PAB.hap1.genes.tsv > AG06213_PAB.hap1.genes.sorted.tsv
sort AG06213_PAB.hap2.genes.tsv > AG06213_PAB.hap2.genes.sorted.tsv
comm -12 AG06213_PAB.hap1.genes.sorted.tsv AG06213_PAB.hap2.genes.sorted.tsv > AG06213_PAB.hap1.hap2.shared_genes.tsv








