/*
=================================================================================================
                        After trimmed Quality-control
=================================================================================================
Here we perform quality control using fastqc
fastqc takes input as reads[0] reads[1] then output logs and html and zip files
*/

// After QC

process AFTERQC {

    tag {"TRIMMED-Data FASTQC on $sample_id"}
    label 'process_low'
    
    publishDir "$projectDir/afterqc" ,  mode:'copy'
    

    input:
    tuple val( sample_id ), path( reads )

    output:
    
    path "afterqc_${sample_id}_logs"

    script:
    """
    mkdir afterqc_${sample_id}_logs
    fastqc -o afterqc_${sample_id}_logs -f fastq -q ${reads}
    """

}

