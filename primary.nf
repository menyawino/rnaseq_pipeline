#!/usr/bin/env nextflow


def currentDir = System.getProperty('user.dir')

// samples_ch = Channel.fromFilePairs("$currentDir/samples/*_R{1,2}_001.fastq.gz", checkIfExists: true)

params.samples = "${launchDir}/samples/*_R{1,2}_001.fastq.gz"

samples_ch = Channel
            .fromPath( params.samples )

process FASTQC {

    tag "running FastQC on ${read.baseName}"

    conda 'workflow/nf/envs/001_QC.yml'

    input:
    path read

    output:
    file 'fastqc_logs'

    script:
    """
    mkdir -p fastqc_logs

    fastqc -o fastqc_logs -f fastq -q ${read}
    """
}

workflow {
    FASTQC(samples_ch)
}

