# RESOURCES_CRC_metacohort_generation

The “RESOURCES_CRC_metacohort_generation” repository contains the scripts that were used to create a metacohort of 1118 colorectal cancer (CRC) patients. This metacohort was generated from public data sources found in Gene Expression Omnibus repository (GEO), specifically, GEO series: GSE14333, GSE143985, GSE17536, GSE17537, GSE33114, GSE37892, GSE38832, GSE39582, GSE92921. All selected datasets were using HG-U133_Plus_2 Array platform from Affymetrics.

### Dataset Information

The CRC_metacohort_July_2024.Rdata file can be found in Zenodo (doi: 10.5281/zenodo.13303050) and contains two data frames enclosing the following information:

1. ‘pheno.data’ data frame includes patients phenotypic information. Rows refer to patients and columns to pheno variables. Disease free survival (DFS) time and event are available for all patients, but other variables might show missing values. 
2. ‘exprs.data’  data frame includes normalized and batch corrected microarray expression data. Rows refer to microarray probes and columns to patients except for the first three colums that indicate corresponding probe id, gene symbol and Jetset score. 

### Folders

The “scripts/” folder contains three R scripts:

1. clinical_data_extraction.R: Script to retrieve patients phenotypic information. Additionally, it filters out patients without DFS time and event info. Duplicated samples are also removed.
2. clinical_data_extraction_misc.R: Script which contains functions that are used in the previous script. It is mandatory to source this scripts to reproduce metacohort generation.
3. complete_metacohort.R: Script to download all the corresponding .CEL files  and performs normalization and batch removal steps.


The “figures/” folder contains a pdf with a summarized description of the CRC metacohort, specifically: (i) workflow of the metacohort generation, (ii) summary table about main phenotypic variables included in the metacohort and (iii) first two principal components before and after applying batch removal.

Please, if you use this CRC metacohort include the Zenodo citation (doi: 10.5281/zenodo.13303050).
