
/*
=================================================================================================
Title : Nextflow workflow on variant analysis
=================================================================================================

Author : Dr. Majeed Jamakhani
=================================================================================================
*/


/*
=================================================================================================
Input Directories
=================================================================================================

*/

params.reads = "s3://mj-nextflow-aws-bucket/data/ggal/*_{1,2}.fq"
params.genome = "s3://mj-nextflow-aws-bucket/genome/transcriptome.fa"


/*
=================================================================================================
Ouput Directories
=================================================================================================

*/

params.trimmed = "s3://mj-nextflow-aws-bucket/results/trimmed"

params.multiqc = "s3://mj-nextflow-aws-bucket/results/multiqc"
params.outdir = "s3://mj-nextflow-aws-bucket/results"

params.qcdir = "s3://mj-nextflow-aws-bucket/results/QC"



/*
=================================================================================================
Channels
=================================================================================================

*/

read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
read_pairs_ch2 = Channel.fromFilePairs(params.reads, checkIfExists: true)
genome_ch = Channel.fromPath(params.genome, checkIfExists: true)



/*
=================================================================================================
Include Modules
=================================================================================================
*/
include {BEFOREQC} from "./modules/beforeqc"
include {SICKLE} from "./modules/sickle"
include {AFTERQC} from "./modules/afterqc"
include {MULTIQC} from "./modules/multiqc"
include {BWAINDEX} from "./modules/bwa_index"


/*
=================================================================================================
                                    Workflow 
=================================================================================================

*/




workflow {

    beforeqc_ch = BEFOREQC (read_pairs_ch)
    
    trimmed_ch = SICKLE (read_pairs_ch)

    afterqc_ch = AFTERQC (trimmed_ch)

    bwa_index_ch = BWAINDEX (genome_ch)


    MULTIQC(beforeqc_ch.mix(afterqc_ch).collect())

}

