#!/usr/bin/env nextflow

params.samplesfile = file("samples.csv")

// Define the path to the sample data file
def currentDir = System.getProperty('user.dir')
params.sampleData = "$currentDir/sample_data.csv"

// Load configuration file
// config = readYaml(file: params.configfile)

workflow {
    main:
        sampleData = prepareSampleData(params.samplesfile)

        def sample_run_ch = Channel
        .fromPath(params.sampleData, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple( row.sample, row.lane, row.read, row.file, row.condition ) }
        // .set { sample_run_ch }

        // sample_mrn = sample_data.unique { it.sample }.collect { it.sample }
        // lane = sample_data.unique { it.lane }.collect { it.lane }
        // read = sample_data.unique { it.read }.collect { it.read }

        fastqc = raw_fastqc()
        // trimmed = trimming(fastqc)
        // post_trim_qc = posttrim_fastqc(trimmed)
        // multiqc = multiqc(fastqc, post_trim_qc)
        // aligned = alignment(trimmed)
        // sorted_bam = sort_bam(aligned)
        // assembled = stringtie_assembly(sorted_bam)

        // if (config['aligner']['tool'] == "hisat2") {
        //     counts = stringtie_count(assembled)
        // } else if (config['aligner']['tool'] == "kallisto") {
        //     counts = kallisto_count(trimmed)
        // } else {
        //     error "Unsupported aligner specified in config.yaml"
        // }

        // if (config['aligner']['tool'] == "hisat2") {
        //     diffexp = ballgown_diffexp(counts)
        // } else if (config['aligner']['tool'] == "kallisto") {
        //     diffexp = sleuth_analysis(counts)
        // }

    // emit:
    //     diffexp
}


process prepareSampleData {
    tag "Preparing sample data"

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


process raw_fastqc {

    tag "Running FastQC on ${sample_mrn}_${lane}_${read}"

    conda "workflow/rules/envs/001_QC.yml"

    cpus 5

    input:
        tuple val(sample), val(lane), val(read), path(file), val(condition) from sample_run_ch

    output:
        path("analysis/001_QC/{sample_mrn}/{sample_mrn}_{lane}_{read}_001_fastqc.html"), 
        path("analysis/001_QC/{sample_mrn}/{sample_mrn}_{lane}_{read}_001_fastqc.zip")
    
    script:
    """
    mkdir -p analysis/001_QC/${sample_mrn}
    fastqc samples/${sample_mrn}_${lane}_${read}_001.fastq.gz -t ${task.cpus} -o analysis/001_QC/${sample_mrn} > logs/001_QC/${sample_mrn}_${lane}_${read}.log 2>&1
    """
}

// process raw_fastqc {
//     conda "workflow/rules/envs/001_QC.yml"

//     input:
//         val(sample_mrn), val(lane), val(read)
    
//     output:
//         tuple val(sample_mrn), val(lane), val(read), path("analysis/001_QC/${sample_mrn}/${sample_mrn}_${lane}_${read}_001_fastqc.html"), path("analysis/001_QC/${sample_mrn}/${sample_mrn}_${lane}_${read}_001_fastqc.zip")
    
//     script:
//     """
//     mkdir -p analysis/001_QC/${sample_mrn}
//     fastqc samples/${sample_mrn}_${lane}_${read}_001.fastq.gz -t ${task.cpus} -o analysis/001_QC/${sample_mrn} > logs/001_QC/${sample_mrn}_${lane}_${read}.log 2>&1
//     """
// }

// process trimming {
//     conda "workflow/rules/envs/001_QC.yml"

//     input:
//         tuple val(sample_mrn), val(lane), val(read), path fastqc_results
    
//     output:
//         tuple val(sample_mrn), val(lane), val(read), path("analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_R1_001_trimmed.fastq.gz"), path("analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_R2_001_trimmed.fastq.gz"), path("analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_fastp.html")
    
//     script:
//     """
//     fastp -i samples/${sample_mrn}_${lane}_R1_001.fastq.gz -I samples/${sample_mrn}_${lane}_R2_001.fastq.gz -o analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_R1_001_trimmed.fastq.gz -O analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_R2_001_trimmed.fastq.gz -h analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_fastp.html --adapter_sequence TACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter_sequence_r2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT > logs/002_trimming/${sample_mrn}_${lane}.log 2>&1
//     """
// }

// process posttrim_fastqc {
//     conda "workflow/rules/envs/001_QC.yml"

//     input:
//         tuple val(sample_mrn), val(lane), val(read), path trimmed_results
    
//     output:
//         tuple val(sample_mrn), val(lane), val(read), path("analysis/003_posttrim_qc/${sample_mrn}/${sample_mrn}_${lane}_${read}_001_trimmed_fastqc.html"), path("analysis/003_posttrim_qc/${sample_mrn}/${sample_mrn}_${lane}_${read}_001_trimmed_fastqc.zip")
    
//     script:
//     """
//     mkdir -p analysis/003_posttrim_qc/${sample_mrn}
//     fastqc analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_${read}_001_trimmed.fastq.gz -t ${task.cpus} -o analysis/003_posttrim_qc/${sample_mrn} > logs/003_posttrim_qc/${sample_mrn}_${lane}_${read}.log 2>&1
//     """
// }

// process multiqc {
//     input:
//         tuple val(sample_mrn), val(lane), val(read), path fastqc_results, path post_trim_qc_results
//     output:
//         path "analysis/003_posttrim_qc/multiqc_raw"
//     script:
//     """
//     multiqc analysis/001_QC/ -o analysis/003_posttrim_qc/multiqc_raw > logs/003_posttrim_qc/multiqc_raw.log 2>&1
//     """
// }



// process multiqc {
//     input:
//     tuple val(sample_mrn), val(lane), val(read), path fastqc_results, path post_trim_qc_results
//     output:
//     path "analysis/003_posttrim_qc/multiqc_raw"
    
//     script:
//     """
//     multiqc analysis/001_QC/ -o analysis/003_posttrim_qc/multiqc_raw > logs/003_posttrim_qc/multiqc_raw.log 2>&1
//     """
// }

// process alignment {
//     input:
//         tuple val(sample_mrn), val(lane), path trimmed_results
//     output:
//         path "analysis/004_alignment/hisat2/${sample_mrn}_${lane}/${sample_mrn}_${lane}.bam"
//     script:
//     """
//     hisat2 --dta --summary analysis/004_alignment/hisat2/${sample_mrn}_${lane}/${sample_mrn}_${lane}.bam.summary -p ${task.cpus} -x ${params.index_hisat2} -1 analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_R1_001_trimmed.fastq.gz -2 analysis/002_trimming/${sample_mrn}/${sample_mrn}_${lane}_R2_001_trimmed.fastq.gz -S analysis/004_alignment/hisat2/${sample_mrn}_${lane}/${sample_mrn}_${lane}.bam > logs/004_alignment/${sample_mrn}_${lane}_alignment.log 2>&1
//     """
// }

// process sort_bam {
//     input:
//         tuple val(sample_mrn), val(lane), path aligned_results
//     output:
//         path "analysis/004_alignment/${sample_mrn}/${sample_mrn}_${lane}_Aligned.sortedByCoord.out.bam"
//     script:
//     """
//     samtools sort -@ ${task.cpus} -o analysis/004_alignment/${sample_mrn}/${sample_mrn}_${lane}_Aligned.sortedByCoord.out.bam ${aligned_results} > logs/004_alignment/${sample_mrn}_${lane}_sort.log 2>&1
//     """
// }

// process stringtie_assembly {
//     input:
//         tuple val(sample_mrn), val(lane), path sorted_bam_results
//     output:
//         path "analysis/005_assembly/${sample_mrn}/${sample_mrn}_${lane}.gtf"
//     script:
//     """
//     stringtie ${sorted_bam_results} -o analysis/005_assembly/${sample_mrn}/${sample_mrn}_${lane}.gtf -p ${task.cpus} 2> logs/005_stringtie/${sample_mrn}_${lane}.log
//     """
// }

// process stringtie_merge {
//     input:
//         path( expand("analysis/005_assembly/{sample}/{sample}_{lane}.gtf", sample: sample_mrn, lane: lane) )
//     output:
//         path "analysis/005_assembly/merged.gtf"
//     script:
//     """
//     stringtie --merge -o analysis/005_assembly/merged.gtf -p ${task.cpus} ${input} 2> logs/005_stringtie/merged.log
//     """
// }

// process stringtie_count {
//     input:
//         tuple val(sample_mrn), val(lane), path sorted_bam_results, path merged_gtf
//     output:
//         path "analysis/006_count/stringtie/${sample_mrn}_${lane}/${sample_mrn}_${lane}.counts"
//     script:
//     """
//     stringtie -e -B -p ${task.cpus} -G ${merged_gtf} -o analysis/006_count/stringtie/${sample_mrn}_${lane}/${sample_mrn}_${lane}.counts ${sorted_bam_results} 2> logs/006_count/stringtie/${sample_mrn}_${lane}.log
//     """
// }

// process kallisto_count {
//     input:
//         tuple val(sample_mrn), val(lane), path trimmed_results
//     output:
//         path "analysis/006_count/kallisto/${sample_mrn}_${lane}"
//     script:
//     """
//     kallisto quant -i ${params.index_kallisto} -o analysis/006_count/kallisto/${sample_mrn}_${lane} -t ${task.cpus} -b 100 ${trimmed_results} > logs/006_count/kallisto/${sample_mrn}_${lane}.log 2>&1
//     """
// }

// process ballgown_diffexp {
//     input:
//         tuple val(sample_mrn), val(lane), path stringtie_counts
//     output:
//         path "results/diffexp/stringtie/${sample_mrn}_${lane}.diffexp.txt"
//     script:
//     """
//     Rscript scripts/ballgown.R --input ${stringtie_counts} --output results/diffexp/stringtie/${sample_mrn}_${lane}.diffexp.txt --log logs/diffexp/${sample_mrn}_${lane}.log
//     """
// }

// process sleuth_analysis {
//     input:
//         tuple val(sample_mrn), val(lane), path kallisto_counts
//     output:
//         path "results/sleuth/differential_expression_results.tsv"
//     script:
//     """
//     Rscript workflow/scripts/sleuth_test.R --output results/sleuth/differential_expression_results.tsv --log logs/007_sleuth_analysis.log
//     """
// }
