# script to find which allele in the atlantic is colsest to the DM allele

cd /scratch/nadjafn/Atlantic_ASE/DM_allele_atlantic_BLAST

# make a BLAST all vs all of the spliced transcripts
gff_atlantic=/scratch/nadjafn/Atlantic_ASE/output/liftoff_syntelogs_v6/updated_v62Atl_liftoff_a0.9s0.9_ALL.gff



conda activate gffread
# extract the spliced transcripts from atlantic
gffread -w atlantic_spliced.cdna.fa -g /scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa $gff_atlantic


# get the spliced transcripts from the DM
gff_DM=/scratch/nadjafn/reference/v6.1/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3
fasta_DM=/scratch/nadjafn/reference/v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa

# extract the spliced transcripts from DM
gffread -w DM_spliced.cdna.fa -g $fasta_DM $gff_DM


conda activate blast
# make a blast database
makeblastdb -in DM_spliced.cdna.fa -dbtype nucl -out DM_spliced.cdna.fa.db

# run the blast
blastn -query atlantic_spliced.cdna.fa -db DM_spliced.cdna.fa.db -outfmt 6 -out atlantic_DM_spliced.cdna.fa.blastn


# add the header 
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > header
cat header atlantic_DM_spliced.cdna.fa.blastn > atlantic_DM_spliced.cdna.fa.blastn.header


# Also for vizualisation align atlantic towards DM spliced transcripts with minimap

conda activate minimap2

minimap2 -ax splice -uf --secondary=yes DM_spliced.cdna.fa atlantic_spliced.cdna.fa > atlantic_2_DM_minimap.sam

# filter reads with MapQ < 30
samtools view -q 30 -bS atlantic_2_DM_minimap.sam > atlantic_2_DM_minimap.bam

# sort the bam file
samtools sort atlantic_2_DM_minimap.bam > atlantic_2_DM_minimap.sorted.bam
# index the bam file
samtools index atlantic_2_DM_minimap.sorted.bam

# remove sam file
rm atlantic_2_DM_minimap.sam


