# arabidopsis_transcriptome_analysis
**PURPOSE**

The primary purpose of this project is to conduct a comprehensive transcriptomic analysis of *Arabidopsis thaliana* under drought stress conditions using RNA-Seq data. The goal is to understand the molecular response of plants to abiotic stress by identifying differentially expressed genes (DEGs) that are activated or suppressed during water deficiency. By doing so, we can pinpoint specific genes, pathways, and biological processes that are involved in drought tolerance. These insights can ultimately support plant biotechnology, breeding programs, and agricultural research aimed at developing drought-resistant crops.


**INTRODUCTION**

Drought is one of the most severe abiotic stresses impacting global agriculture, affecting crop yield, plant development, and overall food security. As climate change intensifies, it becomes increasingly important to understand how plants respond at the molecular level to such environmental challenges.

*Arabidopsis thaliana*, a model organism in plant biology, provides a well-annotated genome and numerous genetic tools, making it ideal for stress-response studies. RNA sequencing (RNA-Seq) is a powerful next-generation sequencing (NGS) technique that allows researchers to profile the entire transcriptome, measure gene expression levels, and detect novel transcripts with high resolution and accuracy.

In this project, we have employed an end-to-end RNA-Seq data analysis pipeline — starting from raw FASTQ files to final DEG identification — to investigate how *Arabidopsis thaliana* modulates gene expression in response to drought stress. The analysis includes quality control, read alignment to the reference genome, transcript assembly, quantification, differential expression analysis, and data visualization.

By comparing gene expression profiles between drought-treated and control samples, we aim to:

- Identify key genes upregulated or downregulated during stress.
- Explore stress-responsive pathways such as ABA signaling, oxidative stress response, and osmoprotection.
- Build a foundation for further functional genomics and systems biology analysis.

**DEINITION**

This project is a hands-on bioinformatics analysis designed to simulate a real-world research pipeline in plant molecular biology. It involves the application of various open-source tools and programming languages (mainly R and command-line tools) to process and interpret RNA-Seq data. Below is a breakdown of the major components:

- **Data Acquisition**: Downloading RNA-Seq datasets from public repositories such as NCBI SRA or ENA for *Arabidopsis thaliana* exposed to control and drought conditions.
- **Quality Control**: Using FastQC to assess the quality of raw reads and determine if trimming or filtering is necessary.
- **Read Alignment**: Aligning reads to the *Arabidopsis thaliana* reference genome using HISAT2

**RECOMMENDED ENVIRONMENT SETUP**

we can install all required tools using Conda:

<pre>bash
conda create -n rnaseq_env python=3.10
conda activate rnaseq_env </pre>

**TOOLS**
| Tool           | Purpose                                 |
| -------------- | --------------------------------------- |
| FastQC         | Raw read quality control                |
| HISAT2         | Spliced alignment to reference genome   |
| SAMtools       | File manipulation and format conversion |
| StringTie      | Transcript assembly and quantification  |
| DESeq2 (R)     | Differential gene expression analysis   |
| MultiQC        | Aggregate and visualize QC reports      |
| IGV            | Visualization of mapped reads           |

**FILE FORMAT USED**
| File Format      | Description                            |
| ---------------- | -------------------------------------- |
| `.fastq.gz`      | Raw sequencing reads                   |
| `.sam` / `.bam`  | Aligned read files                     |
| `.gtf`           | Gene annotation / transcript assembly  |
| `.txt` / `.csv`  | DEG result tables and count matrices   |
| `.html` / `.pdf` | QC and analysis reports                |
| `.r`             | R scripts for DESeq2 and visualization |

**WORKFLOW**
1. Data Acquisition
2. Quality Control
3. Genome Alignment
4. Transcript Assembly and Quantification
5. Differential Gene Expression Analysis
6. Visualization and Interpretation

**STEP 1. DATA ACQUISTION**

**SAMPLE AND REFERENCE GENOME**

This project uses publicly available RNA-Seq datasets and reference genome files for Arabidopsis thaliana to analyze transcriptional responses under drought stress. 

**RNA-Seq Sample Data**

The following paired-end samples were downloaded from the European Nucleotide Archive (ENA): 

-Control Sample: SRR5489633 

-Drought-Stressed Sample: SRR5489636 

<pre> bash
# Control Sample (SRR5489633)
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR548/003/SRR5489633/SRR5489633_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR548/003/SRR5489633/SRR5489633_2.fastq.gz

# Drought Stress Sample (SRR5489636)
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR548/006/SRR5489636/SRR5489636_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR548/006/SRR5489636/SRR5489636_2.fastq.gz </pre>

**Reference Genome and Annotation (TAIR10)**

Reference files were downloaded from Ensembl Plants:

<pre>bash
# Reference Genome FASTA
wget ftp://ftp.ensemblgenomes.org/pub/release-59/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Gene Annotation GTF
wget ftp://ftp.ensemblgenomes.org/pub/release-59/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz

# Unzip downloaded files
gunzip *.gz </pre>

**STEP 2.  QUALITY CONTROL**

Check the quality of the raw FASTQ files using FastQC and summarize all reports with MultiQC.

 Directory: qc_reports/

<pre>bash
mkdir qc_reports && cd qc_reports
fastqc ../data/raw_reads/*.fastq.gz
multiqc . </pre>

Output: Individual .html QC reports and a combined MultiQC report.

**STEP 3. GENOME ALIGNMENT**

Align clean reads to the reference genome using HISAT2 and convert, sort, and index with SAMtools.

Directory: alignment/

##Index the reference genome:
<pre>bash
hisat2-build ../data/ref_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa arabidopsis_index</pre>

##Align reads:
<pre>bash
hisat2 -x arabidopsis_index -1 SRR5489633_1.fastq.gz -2 SRR5489633_2.fastq.gz -S control.sam
hisat2 -x arabidopsis_index -1 SRR5489636_1.fastq.gz -2 SRR5489636_2.fastq.gz -S drought.sam </pre>

##Convert and sort BAM files:
<pre>bash
samtools view -Sb control.sam > control.bam
samtools sort control.bam -o control.sorted.bam
samtools index control.sorted.bam

samtools view -Sb drought.sam > drought.bam
samtools sort drought.bam -o drought.sorted.bam
samtools index drought.sorted.bam </pre>

**STEP4. TRANSCRIPT ASSEMBLY AND QUANTIFICATION**

Use StringTie to assemble transcripts and estimate their expression levels. Then generate a count matrix using the provided script prepDE.py.

Directory: quantification/

<pre>bash
stringtie control.sorted.bam -G ../data/ref_genome/Arabidopsis_thaliana.TAIR10.59.gtf -o control.gtf -e -B
stringtie drought.sorted.bam -G ../data/ref_genome/Arabidopsis_thaliana.TAIR10.59.gtf -o drought.gtf -e -B </pre>

Then run prepDE.py to extract read counts:

<pre>bash
python prepDE.py -i sample_list.txt </pre>

Sample sample_list.txt content:

<pre>bash
control	control/control.gtf
drought	drought/drought.gtf </pre>

**STEP 5. DIFFERENTIAL GENE EXPRESSION ANALYSIS** 

Perform DEG analysis using DESeq2 in R, based on the count matrix generated.

Directory: DEG_analysis/

<pre>
library(DESeq2)

# Load count data and metadata
countData <- read.csv("gene_count_matrix.csv", row.names=1)
colData <- data.frame(
  condition = factor(c("control", "drought")),
  row.names = colnames(countData)
)

dds <- DESeqDataSetFromMatrix(countData, colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

# Order by adjusted p-value
resOrdered <- res[order(res$padj),]

# Export DEGs
write.csv(as.data.frame(resOrdered), "DEGs_results.csv") </pre>

**STEP 6. VISUALIZATION AND INTERPRETATION**  

Visualize expression differences using PCA, Heatmaps, and Volcano plots in R.

Directory: visualizations/


<pre>
library(pheatmap)
library(EnhancedVolcano)

# Heatmap of top 30 genes
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:30]
pheatmap(assay(dds)[select,], cluster_rows=TRUE)

# Volcano Plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj")  </pre>


**CONCLUSION**

This project successfully demonstrated a complete RNA-Seq analysis pipeline—from raw data acquisition to biological interpretation of differentially expressed genes—in Arabidopsis thaliana under drought stress. Using robust tools such as HISAT2, StringTie, and DESeq2, we were able to:
Identify genes that are significantly upregulated or downregulated during drought conditions.
Understand the molecular response of Arabidopsis to abiotic stress.
Highlight potential stress-responsive genes that could be valuable targets for further functional studies or genetic engineering to improve drought tolerance.
This workflow is adaptable and can be extended to study transcriptomic responses in other organisms or under different experimental conditions.

**KEY LEARNINGS** 
During the execution of this project, the following skills and concepts were developed:

-Data Acquisition from public repositories (ENA, Ensembl Plants)

-Quality Control using FastQC and MultiQC

-Genome Indexing and Alignment with HISAT2

-Transcript Assembly & Quantification using StringTie

-Differential Gene Expression Analysis via DESeq2 in R

-Biological Interpretation of transcriptomic data

-Data Visualization using PCA, heatmaps, and volcano plots in R

-Pipeline Documentation and Reproducibility using GitHub

This project also emphasized the importance of open-source tools, automated workflows, and accurate annotation in modern bioinformatics research.





