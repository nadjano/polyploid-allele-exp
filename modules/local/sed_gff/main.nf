process SED_GFF {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(gff)
    val feature

    output:
    tuple val(meta), path("${gff}_${feature}"), emit: gff_generegion

    script:
   
    """
    #!/bin/bash

    # Extract rows with the desired feature and replace the third column with 'exon'
    awk -F"\\t" -v OFS="\\t" -v feature="$feature" '\$3 == feature { \$3 = "exon"; print }' "$gff" > "${gff}_${feature}"
    """
}

