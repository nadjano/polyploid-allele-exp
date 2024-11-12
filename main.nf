
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


params.designfile = "${baseDir}/assets/samples.csv"
params.reference = "${baseDir}/assets/reference.csv"

params.fastq_dir="${baseDir}/fastq_reads"

/*
=================================================================================================
Ouput Directories
=================================================================================================

*/
params.outdir = "${baseDir}/results"


/*
=================================================================================================
Include Modules
=================================================================================================
*/
include {MINIMAP2} from "./modules/minimap2"
include {PAF_TSV} from "./modules/paf_to_tsv"
include {STATS} from "./modules/mapping_stats"
include {PLOT} from "./modules/plot"
include {PLOT_MAPQ_FILTER} from "./modules/plot_mapq_filter"
include {BARPLOT} from "./modules/barplot"
include {PAF_COMPARISON} from "./modules/compare_pafs"
include {BLAST_DB; BLASTN} from "./modules/blastn"
include {MULTIMAPPERS} from "./modules/multimappers"
include {CONCAT} from "./modules/concat_samples_for_experiment"
include {LENGTHBIAS} from "./modules/lengthbias"
include {GFFREAD} from "./modules/nf-core/gffread"
include {MINIMAP2_ALIGN} from "./modules/nf-core/minimap/align"
include {OARFISH} from "./modules/oarfish"
include {BAM2COUNTS} from "./modules/bam2counts"
include {MERGE_COUNTS} from "./modules/mergeCounts"





/*
=================================================================================================
Channels
=================================================================================================

*/



Channel
    .fromPath(params.reference)
    .splitCsv(header: true, sep: ",")
    .map { row -> 
        tuple(id: row.organism,
            row.fasta, 
            row.gtf,
            row.ploidy
        ) 
    }
    .set { reference_ch }

    Channel
    .fromPath(params.designfile)
    .splitCsv(header: true, sep: ",")
    .map { row -> 
        tuple(
            row.organism,
            row.condition, 
            row.replicate,
            row.sample,
            row.fastq_read
        ) 
    }
    .set { samples_ch }
/*
=================================================================================================
                                    Workflow 
=================================================================================================

*/

workflow {

    reference_ch.map {
        id, fasta, gtf, ploidy -> 
        tuple(id, gtf)
    } 
    .set { ch_gtf }

    reference_ch.map {
        id, fasta, gtf, ploidy -> 
        fasta
    } 
    .set { ch_fasta }


    samples_ch
        .map { organism,condition,replicate,sample,fastq_read ->
        tuple(meta = [id: sample, condition: condition, replicate: replicate, organism: organism], fastq_read)
        }
        .set { ch_fastq_reads }


    ch_fastq_reads.view()

    reference_ch.map {
        id, fasta, gtf, ploidy -> fasta
 
    } 
    .set { ch_fasta }

    // ectract gene regions
    gffread_ch = GFFREAD(ch_gtf, ch_fasta)

    // QC for syntelog reference gene lengths
    //ch_gene_lengths = REFERENCE_LENGTHPLOT(gffread_ch.gffread_fasta)

    // map reads to reference
    ch_alignment = MINIMAP2_ALIGN(ch_fastq_reads.combine(gffread_ch.gffread_fasta),
            Channel.value(false),
            Channel.value("bai"),
            Channel.value(false),
            Channel.value(false))

    ch_alignment.bam.view()
    ch_gene_counts = BAM2COUNTS(ch_alignment.bam.combine(gffread_ch.gffread_fasta), Channel.value("false"))
    // collect the gene counts for all samples



    // merged_counts = MERGE_COUNTS(ch_gene_counts.counts.groupTuple(by: 0).view(), Channel.value("true"))

}   

