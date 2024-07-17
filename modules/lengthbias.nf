/*
=================================================================================================
                        ASE counts with length bias 
=================================================================================================
Here we compa
*/

// 

process LENGTHBIAS {

    tag {"Identify length bias between hapl0types $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.length_plot,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val('*pdf')
    script:
    def mm_parameters = mm_parameters.replace(' ', '').replace('-', '')
    """
    module load python
    python3 ${baseDir}/scripts/length_difference.py --paf *competetive*tsv --output_prefix ${experiment}_${mm_parameters}
    """
}