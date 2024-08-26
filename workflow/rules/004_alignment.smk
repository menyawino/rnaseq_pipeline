# A rule to align the trimmed reads to the reference genome using STAR

rule alignment:
    message:
        "Aligning reads to the reference genome"
    input:
        fq1="analysis/002_trimming/{sample}/{sample}_{lane}_R1_001_trimmed.fastq.gz",
        fq2="analysis/002_trimming/{sample}/{sample}_{lane}_R2_001_trimmed.fastq.gz"
    output:
        "analysis/004_alignment/hisat2/{sample}_{lane}/{sample}_{lane}.bam"
    conda:
        "envs/004_alignment.yml"
    threads:
        config["threads"]
    params: 
        path=lambda wildcards: "results/{}".format(wildcards.sample),
        idx=config["aligner"]["index_hisat2"]
    log:
        "logs/004_alignment/{sample}_{lane}_alignment.log"
    benchmark:
        "benchmarks/004_alignment/{sample}/{sample}_{lane}_alignment.txt"
    shell:
        """
        hisat2 \
        --dta \
        --summary {output}.summary \
        -p {threads} \
        -x {params.idx} \
        -1 {input.fq1} \
        -2 {input.fq2} \
        -S {output} \
        > {log} 2>&1
        """


#  A rule to sort the aligned reads with samtools sort

rule sort_bam:
    message: "Sorting the aligned reads"
    input:
        "analysis/004_alignment/hisat2/{sample}_{lane}/{sample}_{lane}.bam"
    output:
        "analysis/004_alignment/{sample}/{sample}_{lane}_Aligned.sortedByCoord.out.bam"
    conda: "envs/004_alignment.yml"
    threads: 8
    params: 
        path=lambda wildcards: "results/{}".format(wildcards.sample)
    log:
        "logs/004_alignment/{sample}_{lane}_sort.log"
    benchmark:
        "benchmarks/004_alignment/{sample}/{sample}_{lane}_sort.txt"
    shell:
        """
        samtools sort \
        -@ {threads} \
        -o {output} \
        {input} \
        > {log} 2>&1
        """
