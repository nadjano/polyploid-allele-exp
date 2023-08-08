
/*
=================================================================================================
Title : Nextflow workflow on variant analysis
=================================================================================================


Author : Dr. Majeed Jamakhani
=================================================================================================
*/

// This is the main workflow


/*
=================================================================================================
Pipeline Input parameters
=================================================================================================

*/

params.reads = "$projectDir/data/*_{1,2}.fq"
params.genome = "$projectDir/genome/genome.fa"
params.known_sites = "$projectDir/known/known_dbsnp138.vcf"
params.adapter_file = "$projectDir/adapter/adapter.fa"

params.trimmed_reads = "$projectDir/trimmed/*_trimmed_{1,2}.fq"

/*
=================================================================================================
Input Channels
=================================================================================================

*/

read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
read_pairs_ch2 = Channel.fromFilePairs(params.reads, checkIfExists: true)
genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
adapter_ch = Channel.fromPath(params.adapter_file, checkIfExists: true)




/*
=================================================================================================
Ouput Directories
=================================================================================================

*/

params.multiqc = "$projectDir/multiqc"
params.outdir = "$projectDir/results"
params.genomedir = "$projectDir/genome"
params.qcdir = "$projectDir/QC"


/*
=================================================================================================
Ouput Channels
=================================================================================================

*/






/*
=================================================================================================
Include Modules
=================================================================================================
*/
include {BEFOREQC} from "./modules/beforeqc"
include {TRIMMOMATIC} from "./modules/trimmomatic"
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

