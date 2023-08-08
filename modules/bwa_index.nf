/*
=================================================================================================
                       bwa index
=================================================================================================
Here we perform genome index using bwa index

*/


params.genomedir = "$projectDir/genome"


process BWAINDEX {

    tag {"bwa index on genome"}

    publishDir params.genomedir, mode:'copy'

    input:
    path genome

    output:
    path ("*")

    script:
    """
    bwa index  $genome 

    """
}

