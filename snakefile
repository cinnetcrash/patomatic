configfile: "config.yaml"

# Function to dynamically get sample names based on the fastq files present in reads_dir
def get_samples():
    import os
    samples = set()
    for filename in os.listdir(config['reads_dir']):
        if filename.endswith("_R1.fastq.gz") or filename.endswith("_R1.fastq"):
            sample_name = filename.split("_R1")[0]
            samples.add(sample_name)
    return list(samples)

# Rule all to specify final targets
rule all:
    input:
        expand("{output_dir}/fastqc/{sample}_R1_fastqc.html", output_dir=config["output_dir"], sample=get_samples()),
        expand("{output_dir}/multiqc/multiqc_report.html", output_dir=config["output_dir"]),
        expand("{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz", output_dir=config["output_dir"], sample=get_samples()),
        expand("{output_dir}/alignment/{sample}_sorted.bam", output_dir=config["output_dir"], sample=get_samples()),
        expand("{output_dir}/qc/{sample}_bam_stats.txt", output_dir=config["output_dir"], sample=get_samples()),
        expand("{output_dir}/vcf/filtered_calls.vcf.gz", output_dir=config["output_dir"]),
        expand("{output_dir}/consensus/consensus_sequence.fasta", output_dir=config["output_dir"]),
        expand("{output_dir}/consensus/quast_report/report.txt", output_dir=config["output_dir"]),
        expand("{output_dir}/prokka/{sample}/prokka_done.txt", output_dir=config["output_dir"], sample=get_samples())

# Example rules adjusted for dynamic sample identification
# FastQC
rule fastqc:
    input:
        lambda wildcards: f"{config['reads_dir']}/{wildcards.sample}_R1.fastq.gz"
    output:
        html="{output_dir}/fastqc/{sample}_R1_fastqc.html",
        zip="{output_dir}/fastqc/{sample}_R1_fastqc.zip"
    shell:
        "fastqc {input} --outdir {output.html.rsplit('/', 1)[0]}"

# MultiQC
rule multiqc:
    input:
        expand("{output_dir}/fastqc/{sample}_R1_fastqc.html", output_dir=config["output_dir"], sample=get_samples())
    output:
        multiqc_report="{output_dir}/multiqc/multiqc_report.html"
    shell:
        "multiqc {input[0].rsplit('/', 2)[0]} -o {output.multiqc_report.rsplit('/', 1)[0]}"

# Trimmomatic for read trimming
rule trim_reads:
    input:
        fwd=lambda wildcards: f"{config['reads_dir']}/{wildcards.sample}_R1.fastq.gz",
        rev=lambda wildcards: f"{config['reads_dir']}/{wildcards.sample}_R2.fastq.gz"
    output:
        trimmed_fwd="{output_dir}/trimmed/{wildcards.sample}_R1_trimmed.fq.gz",
        trimmed_rev="{output_dir}/trimmed/{wildcards.sample}_R2_trimmed.fq.gz"
    shell:
        "trimmomatic PE {input.fwd} {input.rev} {output.trimmed_fwd} /dev/null {output.trimmed_rev} /dev/null ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:{config['trimming']['leading']} TRAILING:{config['trimming']['trailing']} SLIDINGWINDOW:{config['trimming']['slidingwindow']} MINLEN:{config['trimming']['minlen']}"

# BWA for alignment
rule align_reads:
    input:
        fwd="{output_dir}/trimmed/{wildcards.sample}_R1_trimmed.fq.gz",
        rev="{output_dir}/trimmed/{wildcards.sample}_R2_trimmed.fq.gz"
    output:
        bam="{output_dir}/alignment/{wildcards.sample}_sorted.bam"
    shell:
        "bwa mem -t {config['alignment']['threads']} {config['reference']} {input.fwd} {input.rev} | samtools sort -o {output.bam} && samtools index {output.bam}"

# Samtools for generating BAM statistics
rule bam_stats:
    input:
        "{output_dir}/alignment/{wildcards.sample}_sorted.bam"
    output:
        "{output_dir}/qc/{wildcards.sample}_bam_stats.txt"
    shell:
        "samtools stats {input} > {output}"

# BCFtools for variant calling
rule call_variants:
    input:
        bam=expand("{output_dir}/alignment/{sample}_sorted.bam", output_dir=config["output_dir"], sample=get_samples())
    output:
        "{output_dir}/vcf/filtered_calls.vcf.gz"
    shell:
        "bcftools mpileup -Ou -f {config['reference']} {input.bam} | bcftools call -mv -Oz -o {output} && bcftools filter -s LowQual -e '%QUAL<20 or DP<10' {output} -Oz -o {output}"

# BCFtools for generating a consensus sequence
rule generate_consensus:
    input:
        "{output_dir}/vcf/filtered_calls.vcf.gz"
    output:
        "{output_dir}/consensus/consensus_sequence.fasta"
    shell:
        "bcftools consensus -f {config['reference']} {input} > {output}"

# QUAST for quality assessment of the consensus sequence
rule run_quast:
    input:
        "{output_dir}/consensus/consensus_sequence.fasta"
    output:
        report="{output_dir}/consensus/quast_report/report.txt"
    shell:
        "quast {input} -o {output.report.rsplit('/', 1)[0]}"

# Prokka for annotation
rule run_prokka:
    input:
        "{output_dir}/consensus/consensus_sequence.fasta"
    output:
        touch="{output_dir}/prokka/{wildcards.sample}/prokka_done.txt"
    shell:
        "mkdir -p {output.touch.rsplit('/', 1)[0]} && prokka {input} --outdir {output.touch.rsplit('/', 1)[0]} --prefix {wildcards.sample}"
