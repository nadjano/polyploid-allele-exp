/*
=================================================================================================
                        paf files to tsv
=================================================================================================

*/



process PAF_TSV {

    tag {"paf to tsv $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.tsv_out,  mode:'link'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val( strategies ), path( "*.tsv" ), val(mm_parameters)
    script:
    def mm_parameters = mm_parameters.replace(' ', '').replace('-', '')
    """
    module load python
    python3 ${baseDir}/scripts/paf_to_tsv.py ${paf} ${experiment}_${strategies}_${mm_parameters}.tsv
    """
}