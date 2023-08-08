/*
=================================================================================================
                        Multiqc
=================================================================================================
Here we perform quality control using fastqc
fastqc takes input as reads[0] reads[1] then output logs and html and zip files
*/

params.multiqc = "$projectDir/multiqc"

process MULTIQC {

    tag {"Multiqc on $sample_id"}
    publishDir params.multiqc ,  mode:'copy'

    input:
    path ("*")

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

