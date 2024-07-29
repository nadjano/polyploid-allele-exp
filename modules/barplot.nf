/*
=================================================================================================
                        paf files to tsv
=================================================================================================

*/



process BARPLOT {

    tag {"PLot $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.barplot,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val( strategies ), path( "*.pdf" ), val(mm_parameters), path("*tsv")
    script:
    def mm_parameters = mm_parameters.replace(' ', '').replace('-', '')
    """
    module load python
    python3 ${baseDir}/scripts/bar_plot_cats.py \
    --paf ${experiment}_competetive_${mm_parameters}.tsv \
    --output_prefix ${experiment}_${mm_parameters} \
    --min_gene_count 20
    echo "Hello"
    mkdir -p ${experiment}_${mm_parameters}_counts_tsv
    mv *.tsv ${experiment}_${mm_parameters}_counts_tsv/
    
    """
}