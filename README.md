# <img src="https://github.com/leylasrknt/HTS_Analysis/assets/77142451/f22f77ee-7872-42d2-857d-9d237a3decee" alt="HTS Logo" width="80"> HTS_Analysis

A comprehensive collection of scripts and tools for High-Throughput Sequencing (HTS) data analysis. Explore various functionalities, from data preprocessing to differential expression analysis, designed to streamline your omic research projects.

---

# HTS Analysis (BULK RNA SEQ) 

## Run the Script

Execute the script, and it will automatically download the specified SRA file, split it into fastq files, and save them in the specified output directory.

1. **[SRA_download_file.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/SRA_download_file.sh):** This script automates the process of downloading sequencing data from the Sequence Read Archive (SRA), a publicly available repository of high-throughput sequencing data.

2. **[fastqdump.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqdump.sh):** After downloading data from SRA, this script utilizes the SRA Toolkit to convert the downloaded files into FASTQ format. FASTQ is a common file format used to store biological sequences and their corresponding quality scores.

3. **[fastqc_trimmed.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqc_trimmed.sh):** Following data trimming, this script performs quality control using FastQC. FastQC is a tool that analyzes and reports the quality of sequencing data, providing valuable insights into potential issues or artifacts present in the dataset.

4. **[multiqc.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/multiqc.sh):** This script generates a MultiQC report, which compiles and visualizes quality control information from multiple tools used in the analysis. MultiQC simplifies the process of aggregating data from various sources, providing a comprehensive overview of the analysis results.

5. **[trimmer.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/trimmer.sh):** This script trims raw sequencing data, removing low-quality sequences and adapters. Trimming is an essential preprocessing step to enhance the quality of sequencing data, ensuring more accurate downstream analysis.

6. **[hisat2_to_run.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/hisat2_to_run.sh):** Provides instructions and parameters for running HISAT2, a popular tool for aligning sequencing reads to a reference genome or transcriptome. Proper alignment is crucial for understanding where the sequenced fragments originated in the reference genome.

7. **[Picard.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/Picard.sh):** Offers instructions on utilizing the Picard tools, a set of command-line utilities for manipulating high-throughput sequencing data files. Picard tools are commonly used for tasks such as sorting, indexing, and marking duplicates in BAM files.

8. **[bam_to_bami_index.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/bam_to_bami_index.sh):** This script creates index files for BAM (Binary Alignment/Map) files. Indexing BAM files allows for efficient retrieval of specific regions of interest, enabling targeted analysis without the need to process the entire dataset.

9. **[Htseq_to_run.sh](https://github.com/username/repository-name/blob/main/Htseq_to_run.sh):** Provides instructions for running HTSeq, a tool used for counting reads that are mapped to genes. HTSeq is often employed in gene expression analysis, providing gene-level read count information for downstream statistical analysis.

10. **[deseq_timepoint.R](https://github.com/username/repository-name/blob/main/deseq_timepoint.R):** This R script performs differential gene expression analysis using DESeq2 specifically tailored for timepoint data. DESeq2 is a widely used package in R for detecting differentially expressed genes between different conditions or timepoints in RNA-seq data.

11. **[DEseq_heatmap.R](https://github.com/leylasrknt/HTS_Analysis/blob/main/DEseq_heatmap.R):** Generates heatmaps based on the results of DESeq2 analysis. Heatmaps are graphical representations used to visualize gene expression patterns, providing a clear visual representation of differential expression across samples or conditions.

12. **[Upset_Plot.R](https://github.com/leylasrknt/HTS_Analysis/blob/main/Upset_Plot.R):** Creates UpSet plots, a type of data visualization that allows for the analysis of sets, intersections, and their size relationships. UpSet plots are particularly useful for visualizing the overlap and uniqueness of gene sets, providing insights into shared and distinct features across different conditions or datasets.

13. **[GO_clusterProfiler.R](https://github.com/leylasrknt/HTS_Analysis/blob/main/GO_clusterProfiler.R):** This R script performs Gene Ontology (GO) enrichment analysis using the clusterProfiler package. GO enrichment analysis helps identify biological processes, cellular components, and molecular functions that are significantly overrepresented in a set of genes of interest.

Explore, adapt, and utilize these scripts for your HTS data analysis projects. Happy analyzing!

---

## SRA Data Download Script

This Bash script is designed to automate the process of downloading sequencing data from the Sequence Read Archive (SRA). It utilizes the fasterq-dump tool from the NCBI SRA Toolkit to efficiently download SRA files and split them into fastq files.

### How to Use:

Set SRA File:
Replace your_sra_file.sra with the actual name of the SRA file you want to download.

Set Output Directory:
Replace /path/to/your/output_directory with the desired directory where you want to save the downloaded fastq files.

**[SRA_download_file.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/SRA_download_file.sh)**

**[fastqdump.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqdump.sh)**

---

## 🚀FastQC Analysis Scripts

This collection of Bash script is here to make your life easier when it comes to checking the quality of your trimmed FASTQ files using FastQC. The scripts activate the necessary Conda environment, run FastQC on a list of trimmed FASTQ files, and save the results in a specified output directory.

### How to Use:

Replace File Names:
Replace "sample1_trimmed.fastq", "sample2_trimmed.fastq", etc., with the actual filenames of your trimmed FASTQ files.

Set Output Directory:
Replace "/path/to/output_directory" with the desired directory where you want to save the FastQC results.

Run the Script:
Execute the script, and it will automatically run FastQC on the specified trimmed FASTQ files and save the results in the specified output directory.

**[fastqc_trimmed.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/fastqc_trimmed.sh)**

---

## MultiQC Report Generation

Welcome to the MultiQC Report Generation repository! This script simplifies the process of generating comprehensive quality control reports for your bioinformatics analyses.

### What This Script Does:

This Bash script utilizes MultiQC, a powerful tool for aggregating and visualizing quality control information from multiple analysis tools. Specifically, it combines FastQC results from "/results/qc/fastqc" and Trimmomatic logs from "/results/logs/trimmomatic", generating a unified report.

### How to Use:

Specify Input Directories:
Define the paths to your FastQC results and Trimmomatic logs by modifying the fastqc_results and trimmomatic_logs variables, respectively.

Set Output Location:
Customize the multiqc_output variable to specify where you want the MultiQC report to be saved.

Run the Script:
Execute the script, and MultiQC will do the heavy lifting, creating a consolidated report in the specified output directory.

**[multiqc.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/multiqc.sh)**

---

## Trimmomatic Automated Processing Script

Welcome to our Trimmomatic Automated Processing Script repository! This Bash script streamlines the preprocessing of paired-end high-throughput sequencing data using Trimmomatic, a popular tool for adapter removal and quality filtering.

### What This Script Does

This script automates the Trimmomatic process for multiple sets of paired-end sequencing data. It removes adapters, trims low-quality bases, and filters out short sequences, ensuring high-quality data for downstream analysis. The script handles various input file prefixes, making it adaptable to different datasets.

### How to Use

**Specify Input and Output Paths:**
Set the paths for your adapter file, input data directory, output directory, and log files.

**Configure Trimmomatic Parameters:**
Customize the Trimmomatic command according to your data and quality requirements.

**Run the Script:**
Execute the script, and it will automatically process each pair of input files, generating trimmed output files and log reports.

**[trimmer.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/trimmer.sh)**

---

## 🧬HISAT2 Automated Alignment Script

Welcome to our HISAT2 Automated Alignment Script repository! This Bash script facilitates the alignment of paired-end high-throughput sequencing data using HISAT2, a popular tool for mapping sequencing reads to a reference genome.

### What This Script Does

This script automates the alignment process for multiple sets of paired-end sequencing data. It takes trimmed FASTQ files as input, aligns the reads to the specified genome index using HISAT2, converts the output SAM files to sorted BAM files, and indexes them for further analysis. The script is designed to handle various samples, making it adaptable for diverse datasets.

### How to Use

#### Specify Input Paths

Set the paths for your genome index, directory containing trimmed FASTQ files, and the desired output directory.

#### Run the Script

Execute the script, and it will automatically align each pair of trimmed FASTQ files, generating sorted and indexed BAM files.

Feel free to modify the script to match your specific project requirements. If you have any questions or need assistance, don't hesitate to reach out. Enjoy seamless alignment with our automated HISAT2 script!

**[hisat2_to_run.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/hisat2_to_run.sh)**

---

## 🌟 Picard Duplicate Marking Script 

Welcome to our Picard Duplicate Marking Script repository! This Bash script automates the removal of duplicates from your BAM files using Picard, ensuring cleaner and more accurate sequencing data for downstream analysis.

### What This Script Does:

This script processes multiple BAM files, identifies duplicate reads, and marks them for removal. It improves the integrity of your data, preparing it for advanced genomic investigations.

### How to Use:

**Set Paths:** Specify the directories for your input BAM files and where you want to save the Picard results.

**Run the Script:** Execute the script to process each BAM file, removing duplicates and generating marked BAM files along with metrics reports.
**[Picard.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/Picard.sh)**

**GitHub Repository Description:**

🔍 **BAM File Indexing Script**

This Bash script swiftly generates index (BAI) files for your BAM files, ensuring rapid data retrieval and analysis in genomics projects.

### Features:

**Automated Indexing:** Quickly index all BAM files in your specified directory with a single command.

**Customizable Paths:** Easily configure the script by setting your BAM files' directory and the output location for BAI files.

**File Integrity Check:** Ensures data integrity by verifying each file before indexing.

**[bam_to_bami_index.sh](https://github.com/leylasrknt/HTS_Analysis/blob/main/bam_to_bami_index.sh)** 

#### How to Use:

1. **Set BAM File Directory:** Modify `"path/to/your/BAM/files"` to match your directory structure.
2. **Specify Output Directory:** Adjust `"path/to/your/index/directory"` to define the BAI files' storage location.
3. **Run the Script:** Execute the script in your terminal and watch the magic happen!

## HTSeq Counting Script

This Bash script simplifies the process of generating gene counts from your BAM files using HTSeq. Designed for RNA-seq data analysis, this script ensures accurate and efficient counting, providing essential data for downstream analyses.

### Features:
**Automated Counting:** Process multiple BAM files automatically, generating counts for each sample in one go.

**Flexible Input:** Customize input paths for your BAM files and reference annotation file, making it adaptable to various datasets.

**Organized Output:** Counts are saved in separate files, maintaining a clean and organized output directory structure.

### How to Use:
**Set Paths:** Update BAM_DIR, ANNOTATION_FILE, and OUTPUT_DIR to match your file locations.
**Run the Script:** Execute the script, and it will process each BAM file, generating individual count files in the specified output directory.

## 📊 DESeq2 Differential Expression Analysis Script

Welcome to our DESeq2 Differential Expression Analysis Script repository! This R script automates the process of identifying differentially expressed genes from high-throughput sequencing count data, such as RNA-seq counts. It employs the powerful DESeq2 package, streamlining the analysis and providing valuable insights into gene expression variations.

### What This Script Does:

#### Library Loading: 
The script starts by loading essential R libraries, including DESeq2 for statistical analysis and dplyr for data manipulation.

#### Set Your Data: 
Configure the working directory to point to your count data directory, and ensure your metadata is correctly formatted in a CSV file.

#### Data Processing: 
The script scans the directory for count data files with the ".txt" extension, reads, and organizes them. Simultaneously, it reads metadata, preparing the data for analysis.

#### DESeqDataSet Creation: 
The DESeq2 analysis begins by creating a DESeqDataSet, combining count data with metadata. The design formula specifies the variables to analyze, allowing for a targeted investigation.

#### Differential Expression Analysis: 
DESeq2 analysis is performed, identifying genes showing significant expression differences between conditions or time points.

#### Results Export: 
The results are exported as a CSV file ("differential_expression_results.csv"), enabling further exploration and visualization of differentially expressed genes.

#### How to Use:

1. **Modify the working directory path** to match your count data directory.
2. **Ensure your metadata is correctly formatted** in a CSV file and referenced appropriately in the script.
3. **Execute the script** in your R environment. Analyze the generated "differential_expression_results.csv" file to gain valuable insights into gene expression changes.

##   **Optimized Differential Expression Analysis** 

Welcome to our optimized repository for Differential Expression Analysis in R! This code provides efficient scripts for analyzing RNA-seq data using DESeq2. The provided scripts include generating Volcano Plots, PCA Analysis, and Heatmaps. 

**What's Included:**

- **Volcano Plot**: Visualize differential expression results with ease.
- **PCA Analysis**: Explore sample relationships in a high-dimensional space.
- **Heatmap**: Identify gene expression patterns in a visually appealing way.
- **Up/Down-Regulated Genes**: Easily filter genes based on fold change and adjusted p-value thresholds.

**How to Use:**

1. **Clone the Repository**: Download or clone this repository.

2. **Run the Scripts**: Modify parameters in the R scripts and run them to analyze your RNA-seq data.

3. **Analyze Results**: Visualize your data with generated plots and lists.

**[DEseq_heatmap.R](https://github.com/leylasrknt/HTS_Analysis/blob/main/DEseq_heatmap.R)**

# UpSetR Differential Gene Expression Analysis

This R script facilitates the visualization of overlapping and unique gene sets across multiple conditions or time points, helping you gain insights into your RNA-seq data.

### Features:

- **Automated Visualization:** Analyze overlapping gene sets effortlessly with automated UpSet plots.
- **Flexible Input:** Customize input gene lists for different conditions or time points, adapting to your experimental design.
- **Color-Coded Results:** Visualize upregulated and downregulated gene sets with distinct colors, enhancing clarity.
- **Insightful Exploration:** Understand complex gene relationships and identify unique and shared genes across datasets.

### How to Use:

1. **Prepare Gene Lists:**
   - Organize your upregulated and downregulated gene lists into separate files (e.g., `upregulated_genes_4.txt`, `downregulated_genes_4.txt`).
  
2. **Run the Script:**
   - Modify the script to include the paths to your gene list files and customize colors if desired.
   - Execute the script in your R environment to generate UpSet plots for your gene sets.

3. **Explore Your Data:**
   - Analyze the generated UpSet plots to explore overlapping and unique genes among different conditions or time points.
 **[Upset_Plot.R](https://github.com/leylasrknt/HTS_Analysis/blob/main/Upset_Plot.R)**

## Gene Set Enrichment Analysis Script

This R script simplifies the process of performing gene set enrichment analysis using the clusterProfiler package. It enables you to convert gene names to ENSEMBL IDs and identify enriched biological processes (BP) within your gene sets.

Features:
- **Automated Analysis:** Process multiple gene sets automatically, generating enrichment results for each set.
- **Flexible Input:** Customize input gene sets and adjust parameters to match your specific dataset.
- **Organized Output:** Enrichment results are structured and ready for further exploration and interpretation.

How to Use:
1. **Set Your Data:** Replace sample gene lists (e.g., `upregulated_genes_4`, `downregulated_genes_8`) with your gene sets.
2. **Run the Script:** Execute the script in your R environment. Analyze the generated enrichment results to gain insights into biological processes associated with your gene sets.

**[GO_clusterProfiler.R](https://github.com/leylasrknt/HTS_Analysis/blob/main/GO_clusterProfiler.R)**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.

[![License.md: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/leylasrknt/HTS_Analysis/blob/main/LICENSE.md)
