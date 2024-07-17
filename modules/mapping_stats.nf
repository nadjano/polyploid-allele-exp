/*
=================================================================================================
                        concats pafs per experiment
=================================================================================================

*/

process STATS {
    tag {"stats for $sample_id $experiment"}
    label 'process_mem'
  
    publishDir params.stats,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val( strategies ), path( "*.tsv" ), val(mm_parameters)
    script:
    def mm_parameters = mm_parameters.replace(' ', '').replace('-', '')
    """
    module load python
    python3 ${baseDir}/scripts/mapping_stats.py --paf $paf \
    --output_prefix ${experiment}_${sample_id}_${strategies}_${mm_parameters}

    """
}