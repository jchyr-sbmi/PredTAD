# PredTAD

## Getting Started

Dependencies required for PredTAD

Rpackages:
- data.table 
- GenomicRanges
- h2o
- minfi
- plyr
- randomForest
- readr
- ROCR
- ROSE
- rtracklayer
- stringr
- tictoc
- tidyverse

## How to run PredTAD and process data from other cell lines

### Step 1: Data Curation
Create a folder for all cell-specific Hi-C, ChIPseq, histone modification, and methylation data. 
You may use our curated data for non-cell-line specific genomic data, such as distance from centromere (relative chromosomal location), gene density (defined by number of TSS), and number of transcription factor binding sites of 161 transcription factors.
Download, obtain, or generate narrowPeaks bed files for 8 histone modification proteins and 12 transcription factors as listed in List-of-Bed_narrowPeaks files. Create a .txt file with at minimum two columns: file_name and target name. 
Download, obtain, or generate Hi-C data for your cell line. Processed Hi-C data with identified or calculated TAD or insulated regions can also be used. You may use the 3D: Genome Browser to attain TADs for your cell line. 

### Step 2: Data reading and preprocessing
Run data_preprocessing.R to read genomic and epigenomic files. 
TAD boundaries are defined as regions between two adjacent TADs with insulation scores greater than 0.15 (using insulation square analysis at a 40 kb resolution). If TAD regions were used, TAD boundaries are defined as the region half-way between two TADs. Due to variations in replicates, the width of boundaries were extended to 200 kb. 
TF ChIPseq and Histone Modification ChIPseq data are mapped to hg19 with BWA-MEM and narrow peaks files were called with MACS2. They should each contain ten columns: chrom, start, end, name, score, strand, sig, pv, qv, and peak. Average signal values for narrow peaks were used. Values per bin were normalized to 0-1. 
DNA methylation from Illumina Methylation 450K BeadChip array was used in this study. Annotation information for the 450K BeadChip array is provided. Average beta values were used. 

### Step 3: Sample preparation
The entire genome is binned into 10 kb regions. Sample region in the telomere and centromeres were excluded. Average signal values for each bin are normalized from 0 to 1. Ten upstream and ten downstream neighbors’ genomic and epigenomic information were included as feature. 

### Step 4: PredTAD model
Gradient Boosting Machine (GBM) is used as the classifier for PredTAD. The samples of this model are 10 kb genomic regions. The input features are genomic and epigenomic features, such as histone modification, transcription factors, number of TSS, etc. The output is whether or not the region is a TAD boundary or non-TAD boundary. Java’s H2O platform was used to run PredTAD. 

## How to test new samples using our trained model
You can load GBM_model_R_1597829376110_1 for trained model, then test with new data. New data for other cell lines and samples need to be processed in the same manner as the trained data (see above).

## List of files

### Codes 
See Codes folder for PredTAD codes and README_codes.md for description of each file.
- README_codes.md
- 1_gen_pre_info_08132020.R
- 2_sample_preparation.R
- 3_predTAD_model.R
- 4_deg_analysis.R
- 5_feature_distribution_preprocessing.R
- 6_feature_distribution_execution.R

### Data and Results
In the Data folder, you will find a curated list of epigenomic features, including histone modification, methylation, and ATACseq data. 
- List of ChIPseq data
- Trained models
  - GBM_model_R_1591383655864_1
  - GBM_model_R_1591383655864_2512
  - GBM_model_R_1591383655864_7482
  - GBM_model_R_1591383655864_9975
  - GBM_model_R_1591383655864_11903
- Top 5, 10, and 15 features plots

