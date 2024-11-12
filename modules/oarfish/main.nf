process OARFISH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"


    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.meta_info.json")  , emit: info            ,
    optional: false
    tuple val(meta), path("*.quant")  , emit: quant            ,
    optional: false
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def output = "${prefix}"



    """
    oarfish \\
        $args \\
        --output $output \\
        --alignments $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def output = "${prefix}"
    """
    touch $output_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}

