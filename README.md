# BRCA_NA
This repository contains code and pipelines for comparative analysis of **Native American (NA)** and **Caucasian (CA)** breast cancer cohorts across multiple omics layers: DNA mutations, copy number variations (CNV), RNA-seq expression.



## 📂 Repository Structure and Run Order

### 1. DNA Mutation Analysis
Scripts for annotating, converting, filtering, and comparing somatic mutations.

`1. vep_annotation.sh` – Annotate VCF files using VEP  
`2. vcf2maf_batch_convert.sh` – Convert annotated VCF files to MAF format  
`3. filter_and_merge_maf.R` – Filter mutations and merge across cohorts  
`4. fisher_test_and_volcano_plot.R` – Compare mutation frequencies (Fisher’s test) and generate volcano plots  
`UMAP.R` – Visualize mutation profiles using UMAP  



### 2. CNV Analysis
Scripts for CNV calling, filtering, statistical testing, and visualization.

`1. cnv_to_bed_and_annotation.R` – Convert CNV calls to BED and annotate genes  
`2. filter_caucasian_cnv.R` – Filter Caucasian CNV calls  
`3. cnv_statistical_analysis.R` – Compare CNV frequencies across groups  
`4. cnv_volcano_plot.R` – Generate volcano plots for CNV changes  
`5. cnv_results_tables.R` – Summarize CNV results into tables  



### 3. RNA Sequencing Analysis
Scripts for downloading, processing, and analyzing RNA-seq data.

`1. download_rnaseq_caucasian.R` – Download TCGA Caucasian RNA-seq data  
`2. process_na_rnaseq.R` – Process Native American RNA-seq dataset  
`3. combine_rnaseq.R` – Merge NA and CA RNA-seq matrices  
`4. differential_expression.R` – Perform differential expression (DE) analysis  
`5. enrichment_analysis.R` – Pathway enrichment based on DE genes  



### 4. Mutational Signature Analysis
Scripts for extracting and testing mutational signatures.

`1. mutational_signature_pipeline.py` – Run SigProfilerExtractor for de novo mutational signature analysis  
`2. plot&test.py` – Plot exposures and perform statistical tests  



### 5. Multi-Omics Pathway Analysis
Scripts for integrative analysis across DNA, CNV, and RNA.

- `multiomic_kegg_heatmap.R` – Generate KEGG pathway heatmap integrating DNA, CNV, and RNA results  



## ⚙️ Environment and Dependencies

### R
- Main packages (core analysis & visualization):  
  `maftools`, `TCGAbiolinks`, `DESeq2`,  
  `clusterProfiler`, `msigdbr`, `ReactomePA`,  
  `ggplot2`, `ComplexHeatmap`, `pheatmap`  
- Full R session information and package versions:  
  see **`Environment/session_info_R.txt`**

### Python
- Python ≥ 3.10  
- Main packages (mutation signatures & plotting):  
  `SigProfilerExtractor`, `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `biopython`  
- Full Python environment with pinned versions:  
  see **`Environment/requirements_python.txt`**
