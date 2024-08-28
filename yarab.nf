#!/usr/bin/env nextflow

params.samplesfile = file("samples.csv")

// Define the path to the sample data file
def currentDir = System.getProperty('user.dir')
params.sampleData = "$currentDir/sample_data.csv"

workflow {
    main:
        sampleData = prepareSampleData(params.samplesfile)

        Channel
        .fromPath(params.sampleData, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple( row.sample, row.lane, row.read, row.file, row.condition ) }
        .set { sample_run_ch }
        
        sample_run_ch.view()
}


process prepareSampleData {
    conda "workflow/nf/envs/pandas.yml"

    input:
    path samplesfile

    output:
    // path sampleData

    script:
    """
    python3 ${workflow.projectDir}/workflow/nf/scripts/sample_processing.py --csv $samplesfile --dir $projectDir/samples/ 
    """
}