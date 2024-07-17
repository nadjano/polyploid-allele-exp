/*
=================================================================================================
                        Compare paf files
=================================================================================================
Here we compa
*/

// Before QC

process GENE_REGION_ISSUE {

    tag {"Identify problematic gene regions $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.aligned,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( pafs )
    output:    
    tuple  val(experiment), val( sample_id ), 
    script:
    """
    module load python
    python3 ${baseDir}/scripts/
    """
}