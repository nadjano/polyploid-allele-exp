/*
=================================================================================================
                        Long Read Alignment using Minimap2
=================================================================================================
Here we align the long reads to the reference genomes using Minimap2. The alignment is done in two ways:
1. Seperate alignment: In this case, the reads are aligned to each reference genome separately.
2. Competetive alignment: In this case, the reads are aligned to all the reference genomes concatenated together.
    -p 0.8 
                 -N 200

*/              

// After QC

process MINIMAP2 {

    tag {"Minimap on $sample_id"}
    label 'process_minimap2'
    
    publishDir params.aligned,  mode:'link'
    
    input:
    tuple val(experiment),
        val(sample_id), 
        path(reads), 
        val(ploidy), 
        path(gen_ref), 
        path(ref_list), 
        val(strategy),
        val(mm_parameters)

    output:
    tuple  val(experiment), val(sample_id), val(strategy), path("align_${strategy}_${sample_id}*.paf"), val(mm_parameters)

    script:
    def hap_refs_str = ref_list.join(" ")
    def mm_parameters_name = mm_parameters.replace(' ', '').replace('-', '')
    if (strategy == "seperate")
    """
       hap_num=1
        for ref in $hap_refs_str
        do
            module load minimap2
            minimap2 \
                        -t ${task.cpus} \
                        -c \
                        -K 50g \
                        -x splice \
                        -uf \
                        --secondary=yes \
                        ${mm_parameters} \
                        --paf-no-hit \
                        ${params.ref_dir}/\$ref \
                        $reads \
                        > align_${strategy}_${sample_id}_hap\${hap_num}.paf
            hap_num=\$((hap_num+1))
        done
        cat *paf > align_${strategy}_${sample_id}_${mm_parameters_name}.paf
        rm *hap*.paf

    """
    else if (strategy == "competetive") 
        """
        module load minimap2
        minimap2 \
                    -t ${task.cpus} \
                    -c \
                    -K 50g \
                    -x splice \
                    -uf \
                    --secondary=yes \
                    ${mm_parameters} \
                    --paf-no-hit \
                    ${baseDir}/hap_genomes/${gen_ref} \
                    $reads \
                    > align_${strategy}_${sample_id}_${mm_parameters_name}_allhaps.paf
        """
    else
        error "Invalid alignment mode: ${strategy}"
}

//Orangutan_Pf000002_top12000_read_diff.tsv

                       //  -f 0.000002 \
                       // -P 