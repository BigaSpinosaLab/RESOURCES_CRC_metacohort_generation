# RESOURCES_CRC_metacohort_generation

The “RESOURCES_CRC_metacohort_generation” repository contains the scripts that were used to create a metacohort of 1118 colorectal cancer (CRC) patients. This metacohort was generated from public data sources found in Gene Expression Omnibus repository (GEO), specifically, GEO series: GSE14333, GSE143985, GSE17536, GSE17537, GSE33114, GSE38832, GSE39582. All selected datasets were using HG-U133_Plus_2 Array platform from Affymetrics.


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


### References

* *GSE14333*: Jorissen RN, Gibbs P, Christie M, Prakash S et al. Metastasis-Associated Gene Expression Changes Predict Poor Outcomes in Patients with Dukes Stage B and C Colorectal Cancer. Clin Cancer Res 2009 Dec 15;15(24):7642-7651. PMID: 19996206

* *GSE143985*: Shinto E, Yoshida Y, Kajiwara Y, Okamoto K et al. Clinical Significance of a Gene Signature Generated from Tumor Budding Grade in Colon Cancer. Ann Surg Oncol 2020 Oct;27(10):4044-4054. PMID: 32328985

* *GSE17536* and *GSE17537*: Smith JJ, Deane NG, Wu F, Merchant NB et al. Experimentally derived metastasis gene expression profile predicts recurrence and death in patients with colon cancer. Gastroenterology 2010 Mar;138(3):958-68. PMID: 19914252

* *GSE33114*: de Sousa E Melo F, Colak S, Buikhuisen J, Koster J et al. Methylation of cancer-stem-cell-associated Wnt target genes predicts poor prognosis in colorectal cancer patients. Cell Stem Cell 2011 Nov 4;9(5):476-85. PMID: 22056143

* *GSE38832*: Tripathi MK, Deane NG, Zhu J, An H et al. Nuclear factor of activated T-cell activity is associated with metastatic capacity in colon cancer. Cancer Res 2014 Dec 1;74(23):6947-57. PMID: 25320007

* *GSE39582*: Marisa L, de Reyniès A, Duval A, Selves J et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. PLoS Med 2013;10(5):e1001453. PMID: 23700391
