process BAM2STATS {
    tag "$meta.id"

    // replace!!
    conda "/users/nadjafn/.conda/envs/nextflow"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*stats.txt")                  , optional: false, emit: stats
    tuple val(meta), path("*png")                        , optional: false, emit: plot
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${baseDir}/scripts/bam2stats.py \
        --bam  ${bam} \
        --output ${meta.id}
    

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        // python : \$( python --version 2>&1)
    END_VERSIONS
    """

    // stub:
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // def output_file = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    // def bam_index = bam_index_extension ? "touch ${prefix}.bam.${bam_index_extension}" : ""
    // def bam_input = "${reads.extension}".matches('sam|bam|cram')
    // def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    // """
    // touch $output_file
    // ${bam_index}

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     minimap2: \$(minimap2 --version 2>&1)
    // END_VERSIONS
    // """
}









