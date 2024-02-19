#!/bin/bash

# Define variables
REF_GENOME="path/to/reference_genome.fasta"
FORWARD_READS="input_forward.fq.gz"
REVERSE_READS="input_reverse.fq.gz"
ADAPTERS="path/to/adapters.fa"
OUTPUT_DIR="analysis_output"
TRIMMED_DIR="${OUTPUT_DIR}/trimmed_reads"
FASTQC_DIR="${OUTPUT_DIR}/fastqc"
MULTIQC_DIR="${OUTPUT_DIR}/multiqc"
ALIGN_DIR="${OUTPUT_DIR}/alignment"
QC_DIR="${OUTPUT_DIR}/qc"
VCF_DIR="${OUTPUT_DIR}/vcf"
CONS_DIR="${OUTPUT_DIR}/consensus"
PROKKA_DIR="${OUTPUT_DIR}/prokka"

# Create directories
mkdir -p ${OUTPUT_DIR} ${TRIMMED_DIR} ${FASTQC_DIR} ${MULTIQC_DIR} ${ALIGN_DIR} ${QC_DIR} ${VCF_DIR} ${CONS_DIR} ${PROKKA_DIR}

# Step 1: Raw Data Quality Control
echo "Running FastQC..."
fastqc ${FORWARD_READS} ${REVERSE_READS} -o ${FASTQC_DIR}

echo "Running MultiQC..."
multiqc ${FASTQC_DIR} -o ${MULTIQC_DIR}

# Step 2: Trimming and Cleaning
echo "Trimming reads with Trimmomatic..."
trimmomatic PE -phred33 ${FORWARD_READS} ${REVERSE_READS} \
  ${TRIMMED_DIR}/output_forward_paired.fq.gz ${TRIMMED_DIR}/output_forward_unpaired.fq.gz \
  ${TRIMMED_DIR}/output_reverse_paired.fq.gz ${TRIMMED_DIR}/output_reverse_unpaired.fq.gz \
  ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: Alignment to Reference
echo "Aligning reads with BWA..."
bwa mem -t 4 ${REF_GENOME} ${TRIMMED_DIR}/output_forward_paired.fq.gz ${TRIMMED_DIR}/output_reverse_paired.fq.gz > ${ALIGN_DIR}/aligned_reads.sam

# Step 4: Conversion, Sorting, and Indexing of SAM/BAM Files
echo "Processing SAM/BAM files..."
samtools view -bS ${ALIGN_DIR}/aligned_reads.sam > ${ALIGN_DIR}/aligned_reads.bam
samtools sort ${ALIGN_DIR}/aligned_reads.bam -o ${ALIGN_DIR}/sorted_aligned_reads.bam
samtools index ${ALIGN_DIR}/sorted_aligned_reads.bam

# Step 5: BAM File Statistics
echo "Generating BAM file statistics..."
samtools stats ${ALIGN_DIR}/sorted_aligned_reads.bam > ${QC_DIR}/bam_stats.txt
qualimap bamqc -bam ${ALIGN_DIR}/sorted_aligned_reads.bam -outdir ${QC_DIR}/qualimap_report

# Step 6: Variant Calling and Filtration
echo "Calling variants with bcftools..."
bcftools mpileup -Ou -f ${REF_GENOME} ${ALIGN_DIR}/sorted_aligned_reads.bam | bcftools call -mv -Oz -o ${VCF_DIR}/calls.vcf.gz
bcftools filter -s LowQual -e '%QUAL<20 || DP<10' ${VCF_DIR}/calls.vcf.gz > ${VCF_DIR}/filtered_calls.vcf.gz

# Step 7: Consensus Sequence Generation
echo "Generating consensus sequence..."
bcftools consensus -f ${REF_GENOME} ${VCF_DIR}/filtered_calls.vcf.gz > ${CONS_DIR}/consensus_sequence.fasta

# Re-run FastQC and MultiQC on trimmed reads
echo "Running FastQC on trimmed reads..."
fastqc ${TRIMMED_DIR}/*.fq.gz -o ${FASTQC_DIR}/trimmed

echo "Running MultiQC on all FastQC results..."
multiqc ${FASTQC_DIR} -o ${MULTIQC_DIR}

# Step 9: Consensus Sequence Quality Assessment
echo "Assessing consensus sequence quality with QUAST..."
quast.py ${CONS_DIR}/consensus_sequence.fasta -r ${REF_GENOME} -o ${CONS_DIR}/quast_report

# Step 10: Annotation and Functional Analysis (Optional)
echo "Annotating consensus sequence with Prokka..."
prokka ${CONS_DIR}/consensus_sequence.fasta --outdir ${PROKKA_DIR} --prefix annotated

echo "Workflow completed."

