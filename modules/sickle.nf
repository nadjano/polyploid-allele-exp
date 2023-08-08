
/*
    sickle pe -f SRR957824_adapt_R1.fastq -r SRR957824_adapt_R2.fastq \
    -t sanger -o SRR957824_trimmed_R1.fastq -p SRR957824_trimmed_R2.fastq \
    -s /dev/null -q 25
*/
// Trimming


process SICKLE {

    tag {"Sickle Trimming on $sample_id"}
    label 'process_low'
    

    publishDir "$projectDir/trimmed", mode: 'copy'

    input: 
    tuple val(sample_id), path(reads)
    

    output:
    tuple val(sample_id) , path( "*")

    script:
    """
    sickle pe -f ${reads[0]} -r ${reads[1]}  -t sanger -o ${sample_id}_trimmed_1.fq -p ${sample_id}_trimmed_2.fq  -s ${sample_id}_trimmed_combined.fq -q 25
    
    """
}

