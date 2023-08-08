
//  General command to download any files and save to folder is 
//  wget "url-link-to-file" -O outputdir/file.extension     // here -O is capital O

// Data 

/*
 ==============================================

 Reference Download
 url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

 only Mitochondiral region genome 
 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz
 ==============================================
*/
params.genomedir = "$projectDir/genome"

process REFDOWNLOAD{
    publishDir params.genomedir, mode:'copy'

    output: 
        path ("*")

    script:
    """
    mkdir genome

    wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz" -O genome.fa.gz

    """
}



/*
 ==============================================

Known variant sites Data Download

 urls :
 wget  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
 wget  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

 ==============================================
*/

params.knowndir = "$projectDir/known"

process KNOWNDOWNLOAD{
    publishDir params.knowndir, mode:'copy'

    output: 
        path ("*")

    script:
    """
    mkdir known

    wget "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" -O known_dbsnp138.vcf
    
    wget "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx" -O known_dbsnp138.vcf.idx

    """
}




/*
 ==============================================

FASTQ Data Download

FASTQ Data Download

Tumor -Normal data :
urls:

https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz

https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz


 urls :
 
https://zenodo.org/record/1251112/files/raw_child-ds-1.fq
https://zenodo.org/record/1251112/files/raw_child-ds-2.fq

https://zenodo.org/record/1251112/files/raw_mother-ds-1.fq
https://zenodo.org/record/1251112/files/raw_mother-ds-2.fq

 ==============================================
*/

params.data = "$projectDir/data"

process DATADOWNLOAD{
    publishDir params.data, mode:'copy'

    output: 
        path ("*")

    script:
    """
    mkdir data
    

    wget "https://zenodo.org/record/1251112/files/raw_child-ds-1.fq" -O child_1.fq
    
    wget "https://zenodo.org/record/1251112/files/raw_child-ds-2.fq" -O child_2.fq

    wget "https://zenodo.org/record/1251112/files/raw_mother-ds-1.fq" -O mother_1.fq
    
    wget "https://zenodo.org/record/1251112/files/raw_mother-ds-2.fq" -O mother_2.fq



    """
}


/*

Adapter sequence download
url:
https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa
*/

params.adapter = "$projectDir/adapter"

process ADAPTERDOWNLOAD {
    publishDir params.adapter , mode:'copy'

    output: 
        path ("*")

    script:
    """
    wget "https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa" -O adapter.fa

    """
}


/*
 ==============================================

 WorkFLow
 ==============================================
*/


workflow {

    REFDOWNLOAD()
    KNOWNDOWNLOAD()
    DATADOWNLOAD()
    ADAPTERDOWNLOAD()
}

