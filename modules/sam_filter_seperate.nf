/*
=================================================================================================
                        Multiqc
=================================================================================================
Here we perform quality control using fastqc
fastqc takes input as reads[0] reads[1] then output logs and html and zip files
*/



process SAM_FILTER_SEP {

    tag {"Sep $sample_id"}
    publishDir params.aligned ,  mode:'copy'

    input:
    tuple val(sample_id), val(strategy), path(alignment)

    output:
    path "*"

    script:
    def sam_files = alignment.join(" ")
    """
    echo $sam_files
    module load samtools
    # Convert the string back to an array
    IFS=' ' read -r -a array <<< "$sam_files"
    # Select reads based on MAPQ score
    (set -x; samtools view \${array[0]} | \
        cut -f1,5 | sort -k1,1 > mapq_hap1.txt)
    (set -x; samtools view \${array[1]} | \
        cut -f1,5 | sort -k1,1 > mapq_hap2.txt)


    (set -x; join -1 1 -2 1 mapq_hap1.txt mapq_hap2.txt -e 0 > mapq_hap.txt)

    (set -x; awk '\$2 > \$3' mapq_hap.txt | cut -f1 -d' ' > reads_from_hap1.txt)
    echo "Hap1 reads" \$(wc -l reads_from_hap1.txt | cut -f1 -d' ') >> ../hap_aware.log
    (set -x; awk '\$2 < \$3' mapq_hap.txt | cut -f1 -d' ' > reads_from_hap2.txt)
    echo "Hap2 reads" \$(wc -l reads_from_hap2.txt | cut -f1 -d' ') >> ../hap_aware.log
    (set -x; awk '\$2 == \$3' mapq_hap.txt | cut -f1 -d' ' > reads_equal.txt)
    echo "Equal reads" \$(wc -l reads_equal.txt | cut -f1 -d' ') >> ../hap_aware.log

    # (set -x; sort -R reads_equal.txt > temp.txt)
    # (set -x; head -n \$((\$(wc -l reads_equal.txt | cut -f1 -d' ') / 2)) reads_equal.txt >> reads_from_hap1.txt)
    # (set -x; tail -n +\$((\$(wc -l reads_equal.txt | cut -f1 -d' ') / 2 + 1)) reads_equal.txt >> reads_from_hap2.txt)

    #Compile new bam file
    (set -x; samtools view \${array[0]} | grep -wFf reads_from_hap1.txt - > reads_aln_sorted.hap1_select.sam)
    (set -x; samtools view \${array[1]} | grep -wFf reads_from_hap2.txt - > reads_aln_sorted.hap2_select.sam)

    (set -x; samtools view -H \${array[0] > header.txt)
    cat header.txt reads_aln_sorted.hap1_select.sam reads_aln_sorted.hap2_select.sam
    """
}

