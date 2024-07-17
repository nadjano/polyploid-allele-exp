/*
=================================================================================================
                        Long Read Alignment using Minimap2
=================================================================================================
Here we align the long reads to the reference genomes using Minimap2. The alignment is done in two ways:
1. Seperate alignment: In this case, the reads are aligned to each reference genome separately.
2. Competetive alignment: In this case, the reads are aligned to all the reference genomes concatenated together.
*/

// After QC

process BLASTN {

    tag {"BLASTN on $sample_id"}
    label 'process_high'
    
    publishDir params.blast_out,  mode:'copy'
    
    input:
    tuple val(experiment),
        val(sample_id), 
        path(reads), 
        val(ploidy), 
        path(gen_ref), 
        path(ref_list), 
        path(blast_db)

    output:
    tuple  val(experiment), val(sample_id), path("*BLAST.tsv")

    script:
    """
    module load samtools/1.18 ncbi_blast/2.14.1

    # convert fastq to fasta
    awk '{if(NR%4==1) print ">"\$1; else if(NR%4==2) print \$1}' $reads > ${sample_id}.fasta


    ## blast to database
    blastn -db $blast_db/${experiment}_blast \
        -gapopen 5 \
        -gapextend 2 \
        -reward 1 \
        -penalty -1 \
        -num_threads ${task.cpus} \
        -query ${sample_id}.fasta \
        -outfmt "6 qseqid stitle qstart qend qlen sstart send slen length nident bitscore evalue mismatch gaps" \
        -max_target_seqs 200 \
        > ${experiment}_${sample_id}.BLAST.tsv

    echo \$'qseqid\\tstitle\\tqstart\\tqend\\tqlen\\tsstart\\tssend\\tslen\\tlength\\tnident\\tbitscore\\tevalue\\tmismatch\\tgaps' | \
        cat - ${experiment}_${sample_id}.BLAST.tsv > ${experiment}_${sample_id}_head.BLAST.tsv
    rm ${experiment}_${sample_id}.BLAST.tsv
    """
}

process BLAST_DB {

    tag {"Building blast db for $experiment"}
    label 'process_high'
    
    
    input:
    tuple val(experiment), 
        path(gen_ref)

    output:
        tuple val(experiment),
        path("blast_df_${experiment}")

    script:
    """
    module load ncbi_blast/2.14.1
    BLASTDB=blast_df_${experiment}/${experiment}_blast
    mkdir -p blast_df_${experiment}
    makeblastdb -in ${gen_ref} -out \$BLASTDB -dbtype nucl
    """
}
