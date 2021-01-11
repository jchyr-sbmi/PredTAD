## List of files

- 1_gen_pre_info_08132020.R
- 2_sample_preparation.R
- 3_predTAD_model.R
- 4_deg_analysis.R
- 5_feature_distribution_preprocessing.R
- 6_feature_distribution_execution.R

## 1_gen_pre_info_08132020.R 
This file preprocesses epigenomic and genomic data such as Hi-C, ChIPseq, histone modification, methylation data, as well as distance from centromere (relative chromosomal location), gene density (defined by number of TSS), and number of transcription factor binding sites of 161 transcription factors. 

This file also preprocesses Hi-C data for TAD boundary information.

## 2_sample_preparation.R
Next, this file bins the entire genome into 10 kbp regions. Each genomic region will represent one sample. Features for each sample are the preprocessed epigenomic and genomic features generated in the previous step. 

## 3_predTAD_model.R
This file contains the codes for PredTAD. The classifier used is Gradient Boosting (GBM). The samples of this model are 10 kb genomic regions. The input features are genomic and epigenomic features, such as histone modification, transcription factors, number of TSS, etc for the 10 kbp region as well as ten upstream and ten downstream bins. The output is whether or not the region is a TAD boundary or non-TAD boundary. Java’s H2O platform was used to run PredTAD. 

## 4_deg_analysis.R
This file analyzes differentially expressed genes (DEGs) between normal and breast cancer cell lines. Additionally, this file also identifies DEGs that are estrogen-related genes, oncogenes, or survival-related genes. Moreover, using enhancer-promoter information, this file also checks for DEGs with enhancers within or across TAD boundaries. 

## 5_feature_distribution_preprocessing.R
TAD boundaries between normal and cancer cell lines can be shared/common, or unique to one or the other cell line. This file first distinguishes whether a TAD boundary is gained, lost, overlapping, or completely overlapping. Then, epigenomic and genomic features are extracted for each TAD boundary. Next, this file generates training and testing datasets for TAD boundary alteration modeling.

## 6_feature_distribution_execution.R
This file contains the codes to model TAD boundary alterations and to analyze the feature importance. 
