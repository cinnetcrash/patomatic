# PatoMatic
## Bioinformatics Workflow for Illumina Data Analysis

This Snakemake workflow provides a comprehensive pipeline for processing Illumina sequencing data, from quality control through to consensus sequence generation and annotation.

## Workflow Steps

1. **Quality Control**: Assess the quality of raw FASTQ files using FastQC.
2. **MultiQC Report Generation**: Compile FastQC reports into a single MultiQC report for easy visualization.
3. **Read Trimming**: Trim adapters and low-quality bases from reads using Trimmomatic.
4. **Read Alignment**: Align trimmed reads to a reference genome using BWA.
5. **Alignment Sorting and Indexing**: Process-aligned reads with SAMtools to sort and index BAM files.
6. **BAM File Statistics**: Generate statistics for BAM files using SAMtools.
7. **Variant Calling**: Call variants using bcftools.
8. **Consensus Sequence Generation**: Generate a consensus sequence from the variant calls.
9. **Quality Assessment of Consensus Sequence**: Evaluate the consensus sequence quality with QUAST.
10. **Annotation**: Annotate the consensus sequence using Prokka.

## Dependencies

This workflow requires the following tools, managed via a Conda environment:

- FastQC
- MultiQC
- Trimmomatic
- BWA
- SAMtools
- bcftools
- QualiMap (optional in the provided script, replace with your needs)
- QUAST
- Prokka
- Snakemake

## Setup

### Create the Conda Environment

1. Ensure Conda is installed on your system.
2. Create a Conda environment using the provided `environment.yaml` file:

```bash
conda env create -f environment.yaml
```
Activate the Conda environment:
```bash
conda activate bioinformatics_workflow
```

Configuration
Modify the config.yaml file to list your samples and specify the path to the reference genome.
Place your raw FASTQ files in the designated directory as outlined in the config.yaml.
Running the Workflow
With the Conda environment activated and the configuration set, run the workflow using the following command in the directory containing your Snakefile:
```bash
snakemake --cores all
```

Output
The workflow will generate the following outputs in the designated directories:

Quality control reports in fastqc/ and multiqc/
Trimmed reads in trimmed/
Aligned reads and statistics in alignment/ and qc/
Variant calling and consensus sequences in vcf/ and consensus/
Annotation results in prokka/
