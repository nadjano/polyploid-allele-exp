/*
=================================================================================================
                        ASE counts with length bias 
=================================================================================================
Here we compa
*/

// 

process CAT_BARPLOT {

    tag {"Mapping cats $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.cat_plot,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( paf ), val(mm_parameters)
    output:    
    tuple  val(experiment), val( sample_id ), val('*pdf')
    script:
    def mm_parameters = mm_parameters.replace(' ', '').replace('-', '')
    """
    module load python
    python3 ${baseDir}/scripts/bar_plot_cats.py --paf *competetive*tsv --output_prefix ${experiment}_${mm_parameters}
    """
}