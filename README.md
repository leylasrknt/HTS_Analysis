# <img src="https://github.com/leylasrknt/HTS_Analysis/assets/77142451/f22f77ee-7872-42d2-857d-9d237a3decee" alt="HTS Logo" width="80"> HTS_Analysis

A comprehensive collection of scripts and tools for High-Throughput Sequencing (HTS) data analysis. Explore various functionalities, from data preprocessing to differential expression analysis, designed to streamline your omic research projects.




# HTS Analysis (BULK RNA SEQ) 
Run the Script:
Execute the script, and it will automatically download the specified SRA file, split it into fastq files, and save them in the specified output directory.
1. **[SRA_download_file.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/SRA_download_file.sh):** This script automates the process of downloading sequencing data from the Sequence Read Archive (SRA), a publicly available repository of high-throughput sequencing data.

2. **[fastqdump.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqdump.sh):** After downloading data from SRA, this script utilizes the SRA Toolkit to convert the downloaded files into FASTQ format. FASTQ is a common file format used to store biological sequences and their corresponding quality scores.

3. **[fastqc_trimmed.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqc_trimmed.sh):** Following data trimming, this script performs quality control using FastQC. FastQC is a tool that analyzes and reports the quality of sequencing data, providing valuable insights into potential issues or artifacts present in the dataset.

4. **[multiqc.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/multiqc.sh):** This script generates a MultiQC report, which compiles and visualizes quality control information from multiple tools used in the analysis. MultiQC simplifies the process of aggregating data from various sources, providing a comprehensive overview of the analysis results.

5. **[trimmer.sh](https://github.com/username/repository-name/blob/main/trimmer.sh):** This script trims raw sequencing data, removing low-quality sequences and adapters. Trimming is an essential preprocessing step to enhance the quality of sequencing data, ensuring more accurate downstream analysis.

6. **[hisat2_to_run.sh](https://github.com/username/repository-name/blob/main/hisat2_to_run.sh):** Provides instructions and parameters for running HISAT2, a popular tool for aligning sequencing reads to a reference genome or transcriptome. Proper alignment is crucial for understanding where the sequenced fragments originated in the reference genome.

7. **[Picard.sh](https://github.com/username/repository-name/blob/main/Picard.sh):** Offers instructions on utilizing the Picard tools, a set of command-line utilities for manipulating high-throughput sequencing data files. Picard tools are commonly used for tasks such as sorting, indexing, and marking duplicates in BAM files.

8. **[bam_to_bami_index.sh](https://github.com/username/repository-name/blob/main/bam_to_bami_index.sh):** This script creates index files for BAM (Binary Alignment/Map) files. Indexing BAM files allows for efficient retrieval of specific regions of interest, enabling targeted analysis without the need to process the entire dataset.

9. **[Htseq_to_run.sh](https://github.com/username/repository-name/blob/main/Htseq_to_run.sh):** Provides instructions for running HTSeq, a tool used for counting reads that are mapped to genes. HTSeq is often employed in gene expression analysis, providing gene-level read count information for downstream statistical analysis.

10. **[deseq_timepoint.R](https://github.com/username/repository-name/blob/main/deseq_timepoint.R):** This R script performs differential gene expression analysis using DESeq2 specifically tailored for timepoint data. DESeq2 is a widely used package in R for detecting differentially expressed genes between different conditions or timepoints in RNA-seq data.

11. **[DEseq_heatmap.R](https://github.com/username/repository-name/blob/main/DEseq_heatmap.R):** Generates heatmaps based on the results of DESeq2 analysis. Heatmaps are graphical representations used to visualize gene expression patterns, providing a clear visual representation of differential expression across samples or conditions.

12. **[Upset_Plot.R](https://github.com/username/repository-name/blob/main/Upset_Plot.R):** Creates UpSet plots, a type of data visualization that allows for the analysis of sets, intersections, and their size relationships. UpSet plots are particularly useful for visualizing the overlap and uniqueness of gene sets, providing insights into shared and distinct features across different conditions or datasets.

13. **[GO_clusterProfiler.R](https://github.com/username/repository-name/blob/main/GO_clusterProfiler.R):** This R script performs Gene Ontology (GO) enrichment analysis using the clusterProfiler package. GO enrichment analysis helps identify biological processes, cellular components, and molecular functions that are significantly overrepresented in a set of genes of interest.

Explore, adapt, and utilize these scripts for your HTS data analysis projects. Happy analyzing!

## SRA Data Download Script
This Bash script is designed to automate the process of downloading sequencing data from the Sequence Read Archive (SRA). It utilizes the fasterq-dump tool from the NCBI SRA Toolkit to efficiently download SRA files and split them into fastq files.

### How to Use:
Set SRA File:
Replace your_sra_file.sra with the actual name of the SRA file you want to download.

### Set Output Directory:
Replace /path/to/your/output_directory with the desired directory where you want to save the downloaded fastq files.

**[SRA_download_file.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/SRA_download_file.sh)** 

**[fastqdump.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqdump.sh)** 

## 🚀FastQC Analysis Scripts
 
This collection of Bash script is here to make your life easier when it comes to checking the quality of your trimmed FASTQ files using FastQC. The scripts activate the necessary Conda environment, run FastQC on a list of trimmed FASTQ files, and save the results in a specified output directory.

### How to Use:
Replace File Names:
Replace "sample1_trimmed.fastq", "sample2_trimmed.fastq", etc., with the actual filenames of your trimmed FASTQ files.

### Set Output Directory:
Replace "/path/to/output_directory" with the desired directory where you want to save the FastQC results.

### Run the Script:
Execute the script, and it will automatically run FastQC on the specified trimmed FASTQ files and save the results in the specified output directory.


**[fastqc_trimmed.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqc_trimmed.sh)** 

## MultiQC Report Generation
Welcome to the MultiQC Report Generation repository! This script simplifies the process of generating comprehensive quality control reports for your bioinformatics analyses.

### What This Script Does:
This Bash script utilizes MultiQC, a powerful tool for aggregating and visualizing quality control information from multiple analysis tools. Specifically, it combines FastQC results from "/results/qc/fastqc" and Trimmomatic logs from "/results/logs/trimmomatic", generating a unified report.

### How to Use:
Specify Input Directories:
Define the paths to your FastQC results and Trimmomatic logs by modifying the fastqc_results and trimmomatic_logs variables, respectively.

### Set Output Location:
Customize the multiqc_output variable to specify where you want the MultiQC report to be saved.

### Run the Script:
Execute the script, and MultiQC will do the heavy lifting, creating a consolidated report in the specified output directory.
**[multiqc.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/multiqc.sh)**
