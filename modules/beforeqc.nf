/*
=================================================================================================
                        Quality-control
=================================================================================================
Here we perform quality control using fastqc
fastqc takes input as reads[0] reads[1] then output logs and html and zip files
*/

// Before QC

process BEFOREQC {

    tag {"Raw-Data FASTQC on $sample_id"}
    label 'process_low'
    
    publishDir "$projectDir/beforeqc" ,  mode:'copy'
    

    input:
    tuple val( sample_id ), path( reads )

    output:
    
    path "beforeqc_${sample_id}_logs"

    script:
    """
    mkdir beforeqc_${sample_id}_logs
    fastqc -o beforeqc_${sample_id}_logs -f fastq -q ${reads}
    """

}

