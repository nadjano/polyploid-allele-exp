/*
=================================================================================================
                        concats pafs per experiment
=================================================================================================

*/

process CONCAT {

    tag {"concating $experiment"}
    label 'process_low'
  
    publishDir params.paf_compare,  mode:'link'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val( strategies ), path( "*.paf" ), val(mm_parameters)
    script:
    def mm_parameters_name = mm_parameters.replace(' ', '').replace('-', '')
    """
    cat *paf > ${experiment}_${strategies}_${mm_parameters_name}.paf
    """
}