#!/usr/bin/env nextflow

process runPythonFunction {
    input:
    path csv_file

    output:
    file 'sample_data.csv' into vals

    script:
    """
    python3 get_sample_data.py --csv $csv_file > sample_data.csv
    """
}

workflow {
    csv_file = file('path/to/samples.csv')

    runPythonFunction(csv_file)

    vals.view { file ->
        println "The CSV file created by Python is located at: ${file}"
    }
}
