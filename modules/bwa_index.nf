/*
=================================================================================================
                       bwa index
=================================================================================================
Here we perform genome index using bwa index

*/





process BWAINDEX {

    tag {"bwa index on genome"}

    publishDir params.genome, mode:'copy'

    input:
    path genome

    output:
    path ("*")

    script:
    """
    bwa index  $genome 

    """
}

