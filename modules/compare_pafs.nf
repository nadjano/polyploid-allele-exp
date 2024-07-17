/*
=================================================================================================
                        Compare paf files
=================================================================================================
Here we compa
*/

// Before QC

process PAF_COMPARISON {

    tag {"Compare PAFs for $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.paf_compare,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ), val( strategies ), path( pafs )
    output:    
    tuple  val(experiment), val( sample_id ), path( "*comparison.txt" ), path( "*stats.txt" )
    script:
    """
    module load python
    python3 ${baseDir}/scripts/compare_pafs.py --paf_seperate1 align_seperate_${sample_id}_hap1.paf.tsv  --paf_seperate2 align_seperate_${sample_id}_hap2.paf.tsv  --paf_competetive *competetive*${sample_id}_allhaps.paf.tsv --summary  ${experiment}_${sample_id}_par_comparison.txt --stats ${experiment}_${sample_id}_par_stats.txt
    """
}