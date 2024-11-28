process ALLELEFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "/users/nadjafn/.conda/envs/allelefinder" //conda "${moduleDir}/environment.yml"
    container 'yourdockerhubuser/allelefinder:latest' // Replace with the actual container if available

    input:
        tuple val(meta), path(fasta), path(gff), path(cds)


    output:
        tuple val(meta), path("allelefinder_${identity}ident/allele.adjusted.txt"), emit: construct_output

    script:
        def num_threads = task.cpus ?: 4
        def identity    = task.ext.identity    ?: 90
        def out_dir     = "allelefinder_${identity}ident"

        """
        #!/bin/bash
  
        #echo "export PATH=${params.allelefinder_path}/:\$PATH" >> ~/.bash_profile
        #echo "export PATH=${params.mcscanx_path}/:\$PATH" >> ~/.bash_profile
        source ~/.bash_profile

        # if 'gene' is not present in the 3rd colum gff file, replace 'mRNA' with 'gene'
        sed 's/mRNA/gene/' $gff > ${gff}_gene
        #sed 's/transcript/gene/' $gff > ${gff}_gene

        # only get the lines with chr in the first column
        awk '\$1 ~ /chr/ { print }' ${gff}_gene > ${gff}_chr

        awk '\$1 ~ /_1/ { gsub(/_1/, "A", \$1); print }' ${gff}_chr > ${gff}_chrom
        awk '\$1 ~ /_2/ { gsub(/_2/, "B", \$1); print }' ${gff}_chr >> ${gff}_chrom
        awk '\$1 ~ /_3/ { gsub(/_3/, "C", \$1); print }' ${gff}_chr >> ${gff}_chrom
        awk '\$1 ~ /_4/ { gsub(/_4/, "D", \$1); print }' ${gff}_chr >> ${gff}_chrom


        ${params.allelefinder_path}/allelefinder.py prepare \\
            --in_gff3 ${gff}_chrom  \\
            --in_cds $cds \\
            --out_cds ${cds}_deduplicated \\
            --out_gff3 ${gff}_deduplicated \\


        touch ref.fasta
        touch ref.gff
        touch ref.cds

        echo "starting construction"

       ${params.allelefinder_path}/allelefinder.py construct \\
            --ref ref.fasta \\
            -d ref.gff \\
            -f ref.cds  \\
            -c $cds \\
            -g ${gff}_chrom \\
            -n 4 \\
            -t $num_threads \\
            -w $out_dir \\
            -i $identity

        ${params.allelefinder_path}/allelefinder.py stat \\
            -i $out_dir/allele.adjusted.txt \\
            -g $gff \\
            -o $out_dir

        """
}

