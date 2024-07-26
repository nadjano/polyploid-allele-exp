/*
=================================================================================================
                        paf files to tsv
=================================================================================================

*/



process PLOT {

    tag {"PLot $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.plot,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val( strategies ), path( "*.pdf" ), path( "${experiment}_*_tsv" ), val(mm_parameters)
    script:
    def mm_parameters = mm_parameters.replace(' ', '').replace('-', '')
    """
    module load python
    python3 ${baseDir}/scripts/plot_mapping_strategy_comparison.py \
    --paf_seperate ${experiment}_seperate_${mm_parameters}.tsv \
    --paf_competetive ${experiment}_competetive_${mm_parameters}.tsv \
    --output_prefix ${experiment}_${mm_parameters}

    echo "Hello you!"
    mkdir -p ${experiment}_${mm_parameters}_tsv 
    mv *.tsv ${experiment}_${mm_parameters}_tsv 
    """
}