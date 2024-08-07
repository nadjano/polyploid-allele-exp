
/*
=================================================================================================
Title : Nextflow workflow on mapping strategy comparison for ASE analysis 
=================================================================================================

Author : Nadja Nolte
=================================================================================================
*/


/*
=================================================================================================
Input Directories
=================================================================================================

*/

params.projectdir="/blue/mcintyre/share/potato_ASE/" 
params.designfile = "${params.projectdir}/design_file/nf_mapping_comparison_samples_test.csv"
params.ref_dir="${baseDir}/hap_genomes"
params.fastq_dir="${baseDir}/fastq_reads"

/*
=================================================================================================
Ouput Directories
=================================================================================================

*/

params.aligned = "${baseDir}/out_ms/aligned"
params.paf_compare = "${baseDir}/out_ms/paf_compare"
params.blast_out = "${baseDir}/out_ms/blast_out"
params.multimappers = "${baseDir}/out_ms/mapping_comparison_mm_blast"
params.plot = "${baseDir}/out_ms/scatter_plots"
params.barplot = "${baseDir}/out_ms/bar_plots"
params.stats = "${baseDir}/out_ms/mapping_stats"
params.length_plot = "${baseDir}/out_ms/length_bias"
params.tsv_out = "${baseDir}/out_ms/paf_tsv_files"

/*
=================================================================================================
Channels
=================================================================================================

*/

Channel
    .fromPath(params.designfile)
    .splitCsv(header: true, sep: ",")
    .map { row -> 
        tuple(
            row.experiment,
            row.sample, 
            "${params.fastq_dir}/${row.fastq_read}", 
            row.ploidy, 
            "${params.ref_dir}/${row.gen_ref}", 
            [
                "${params.ref_dir}/${row.hap1_ref}", 
                "${params.ref_dir}/${row.hap2_ref}", 
                "${params.ref_dir}/${row.hap3_ref}", 
                "${params.ref_dir}/${row.hap4_ref}",
                "${params.ref_dir}/${row.hap5_ref}"
            ]
        ) 
    }
    .set { samples_ch }

// Drop the null paths from the samples channel in the hapref list
samples_ch = samples_ch.map { experiment, sample, read, ploidy, gen_ref, hap_refs ->
    tuple(experiment, sample, read, ploidy, gen_ref, hap_refs.findAll { it != null && !it.endsWith("/null") })
}

Channel.fromPath(params.designfile)
    .splitCsv(header: true, sep: ",")
    .map { row -> tuple(Integer.parseInt(row.ploidy)) }
    .set { ploidy_ch }

// Make two channels with mode 'seperate' and 'competetive' for the two different mapping strategies
Channel.from("seperate", "competetive")
    .set { mapping_strategies_ch }

// combine the two channels to get all possible combinations of samples and mapping strategies
samples_ch.combine(mapping_strategies_ch)
    .set { samples_mapping_ch }

// group the samples by experiment
samples_ch
    .groupTuple(by: [0,4])
    .map { experiment, sample, read, ploidy, gen_ref, hap_refs -> tuple(experiment, gen_ref)}
    .set { samples_grouped_ch }



/*
=================================================================================================
Include Modules
=================================================================================================
*/
include {MINIMAP2} from "./modules/minimap2"
include {PAF_TSV} from "./modules/paf_to_tsv"
include {STATS} from "./modules/mapping_stats"
include {PLOT} from "./modules/plot"
include {BARPLOT} from "./modules/barplot"
include {PAF_COMPARISON} from "./modules/compare_pafs"
include {BLAST_DB; BLASTN} from "./modules/blastn"
include {MULTIMAPPERS} from "./modules/multimappers"
include {CONCAT} from "./modules/concat_samples_for_experiment"
include {LENGTHBIAS} from "./modules/lengthbias"


/*
=================================================================================================
                                    Workflow 
=================================================================================================

*/

workflow {
   // mm_paramers = Channel.from("-N 200", "-P")
   mm_paramers = Channel.from("-N 200", "-P", "-P -f 0.000002")
   //mm_paramers = Channel.from("-N 200")
   minimap_out = MINIMAP2(samples_mapping_ch.combine(mm_paramers))
   minimap_out.transpose().view()
   // get mapping stats
   stats_out = STATS(minimap_out)
   stats_out.view()
   // group minimap by experiment and strategy 
   minimap_out
    .groupTuple(by: [0,2,4])
    .set { minimap_out_exp }
   minimap_out_exp.view()
   // concat paf files for each experiment and strategy
   minimap_exp_strat = CONCAT(minimap_out_exp)
   minimap_exp_strat.view()
   
    // convert paf to tsv
   paf_to_tsv = PAF_TSV(minimap_exp_strat)
   paf_to_tsv.view()
   
   paf_to_tsv
    .groupTuple(by: [0,4])
    .set { minimap_out_grouped }
    minimap_out_grouped.view()
    // plot mapping strategy comparison
    PLOT(minimap_out_grouped)
    // plot barplot for ASE count distribution
    BARPLOT(minimap_out_grouped)

    // Filter channel for experiment == Orangutan and Atlantic
    minimap_out_grouped_length = minimap_out_grouped.join(Channel.of("Orangutan", "Atlantic", "Atlantic_withS"))

    LENGTHBIAS(minimap_out_grouped_length)

//     minimap_out_grouped_RIL.view()

//     PAF_COMPARISON(minimap_out_grouped_RIL)

//     minimap_out_grouped_ATLANTIC = minimap_out_grouped.join(Channel.of(["Atlantic"]))
//     minimap_out_grouped_ATLANTIC.view()
    
    // for atlantic use seperate alignments to identify problems with gene regions

}   