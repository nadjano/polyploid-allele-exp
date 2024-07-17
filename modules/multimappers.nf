/*
=================================================================================================
                        Compare paf files
=================================================================================================
Here we compa
*/

// Before QC

process MULTIMAPPERS {

    tag {"Get mapping locations $experiment $sample_id"}
    label 'process_mem'
  
    publishDir params.multimappers,  mode:'copy'
    
    input:
    tuple  val(experiment), val( sample_id ),val(strategy), path( pafs ), path(blast)
    output:    
    tuple val(experiment), val( sample_id ), val(strategy), path("*txt"), path("*stats"), path("*mapping_comparion*"), path('*_agreement')
    script:
    """
    module load python
    # Compare mapping location for mm2 and blast based on AS score
    python3 ${baseDir}/scripts/get_multimappers_from_paf.py --experiment ${experiment} --paf_seperate *sep*${sample_id}.paf --paf_competetive *competetive*${sample_id}.paf --blast *head.BLAST.tsv --summary ${experiment}_${sample_id}_mapping_agreement_summary_AS.txt --mapping ${experiment}_${sample_id}_mapping_comparion_AS.tsv
    # Do the same for NM score
    #python3 ${baseDir}/scripts/get_multimappers_from_paf_NM.py --paf_seperate *sep*paf --paf_competetive *competetive*paf --blast *head.BLAST.tsv --summary ${experiment}_${sample_id}_mapping_agreement_summary_NM.txt --mapping ${experiment}_${sample_id}_mapping_comparion_NM.tsv
    """
}