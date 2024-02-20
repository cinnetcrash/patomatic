configfile: "config.yaml"

SAMPLES, = glob_wildcards(config["reads_dir"] + "/{sample}_R1.fastq.gz")

rule all:
    input:
        fastqc_reports = expand("fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        multiqc_report = "multiqc/multiqc_report.html",
        trimmed_reads = expand("trimmed/{sample}_R1_trimmed.fq.gz", sample=SAMPLES),
        aligned_bam = expand("alignment/{sample}_sorted.bam", sample=SAMPLES),
        bam_stats = expand("qc/{sample}_bam_stats.txt", sample=SAMPLES),
        vcf_file = "vcf/filtered_calls.vcf.gz",
        consensus_sequence = "consensus/consensus_sequence.fasta",
        quast_report = "consensus/quast_report/report.txt",
        prokka_annotation = expand("prokka/{sample}/prokka_done.txt", sample=SAMPLES)

rule fastqc:
    input:
        lambda wildcards: f"{config['reads_dir']}/{wildcards.sample}_R1.fastq.gz"
    output:
        html="fastqc/{sample}_R1_fastqc.html",
        zip="fastqc/{sample}_R1_fastqc.zip"
    conda:
        "envs/environment.yaml" 
    container:
        "staphb/fastqc:latest"
    shell:
        "fastqc {input} --outdir fastqc/"

rule multiqc:
    input:
        expand("fastqc/{sample}_R1_fastqc.html", sample=SAMPLES)
    output:
        "multiqc/multiqc_report.html"
    container:
        "staphb/multiqc:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for MultiQC
    shell:
        "multiqc fastqc/ -o multiqc/"

rule trim_reads:
    input:
        fwd=lambda wildcards: f"{config['reads_dir']}/{wildcards.sample}_R1.fastq.gz",
        rev=lambda wildcards: f"{config['reads_dir']}/{wildcards.sample}_R2.fastq.gz"
    output:
        trimmed_fwd="trimmed/{sample}_R1_trimmed.fq.gz",
        trimmed_rev="trimmed/{sample}_R2_trimmed.fq.gz"
    container:
        "staphb/trimmomatic:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for Trimmomatic
    shell:
        "trimmomatic PE {input.fwd} {input.rev} {output.trimmed_fwd} /dev/null {output.trimmed_rev} /dev/null ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule align_reads:
    input:
        fwd="trimmed/{sample}_R1_trimmed.fq.gz",
        rev="trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        bam="alignment/{sample}_sorted.bam"
    container:
        "staphb/bwa:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for BWA
    shell:
        """
        bwa mem -t 4 {config['reference']} {input.fwd} {input.rev} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule bam_stats:
    input:
        "alignment/{sample}_sorted.bam"
    output:
        "qc/{sample}_bam_stats.txt"
    container:
        "staphb/samtools:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for Samtools
    shell:
        "samtools stats {input} > {output}"

rule call_variants:
    input:
        bam=expand("alignment/{sample}_sorted.bam", sample=SAMPLES)
    output:
        vcf="vcf/filtered_calls.vcf.gz"
    container:
        "staphb/bcftools:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for bcftools
    shell:
        """
        bcftools mpileup -Ou -f {config['reference']} {input.bam} | bcftools call -mv -Oz -o vcf/calls.vcf.gz
        bcftools filter -s LowQual -e '%QUAL<20 or DP<10' vcf/calls.vcf.gz > {output.vcf}
        """

rule generate_consensus:
    input:
        vcf="vcf/filtered_calls.vcf.gz"
    output:
        fasta="consensus/consensus_sequence.fasta"
    container:
        "staphb/bcftools:latest"
    conda:
        "envs/environment.yaml"  # Use the same Conda environment as for calling variants
    shell:
        "bcftools consensus -f {config['reference']} {input.vcf} > {output.fasta}"

rule run_quast:
    input:
        fasta="consensus/consensus_sequence.fasta"
    output:
        report="consensus/quast_report/report.txt"
    container:
        "staphb/quast:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for QUAST
    shell:
        "quast.py {input.fasta} -o consensus/quast_report"

rule run_prokka:
    input:
        fasta="consensus/consensus_sequence.fasta"
    output:
        touch="prokka/{sample}/prokka_done.txt"
    container:
        "staphb/prokka:latest"
    conda:
        "envs/environment.yaml"  # Path to a specific Conda environment file for Prokka
    shell:
        """
        mkdir -p prokka/{wildcards.sample}
        prokka {input.fasta} --outdir prokka/{wildcards.sample} --prefix {wildcards.sample}
        touch {output.touch}
        """
