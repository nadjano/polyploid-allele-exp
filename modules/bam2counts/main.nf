process BAM2COUNTS {
    tag "$meta.id"

    // replace!!
    conda "/users/nadjafn/.conda/envs/nextflow"

    input:
    tuple val(meta), path(bam), val(meta2), path(fasta)
    val filter

    output:
    tuple val(meta), path("*.tsv")                       , optional: false, emit: counts
    tuple val(meta), path("filter/*.json")               , optional: true, emit: json
    tuple val(meta), path("*.counts_stats.txt")   , optional: true, emit: stats
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bam_filtered = "${meta.id}_${meta.condition}.filtered.bam"
    def output_name = "${meta.id}_${meta.condition}.counts.tsv"
    def filter = filter ? "--filter" : ""
    def min_read_length = 300
    """
    # remove reads with read length < min_read_length
    samtools view -m ${min_read_length} $bam -b  > ${bam_filtered}

    python ${baseDir}/scripts/prefilter.py \
        --aln ${bam_filtered} \
        --target ${fasta} \
        --out_dir filter

    python ${baseDir}/scripts/scores_to_genecounts.py \
    --scores "filter/scores.tsv" \
    --sample ${meta.id} \
    --condition ${meta.condition} \
    --output ${output_name}
    
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



