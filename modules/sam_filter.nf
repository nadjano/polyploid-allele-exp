/*
=================================================================================================
                       bwa index
=================================================================================================
Here we perform genome index using bwa index

*/


process SAM_FILTER {

    tag {"Filter $sample_id $strategy"}
    label 'process_high'

    publishDir params.aligned, mode:'copy'

    input:
    tuple val(sample_id), val(strategy), path(alignment)

    output:
    tuple val(sample_id), val(strategy), path("*filtered.sam")

    script:
    def sam_files = alignment.join(" ")
    """
    module load samtools
    # Iterate over sam files in directory
    for file in $sam_files; do
        samtools view -h -F 260 -q 5 \$file > \${file}.filtered.sam
    done
    """
}

