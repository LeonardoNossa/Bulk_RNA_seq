
# RNA Bulk-seq Analysis in R

## Description

The goal of this project is to identify differentially expressed genes (DEGs) between multiple conditions in bulk RNA-Seq data. The data were collected from [source], and the analysis focuses on the following steps:

- Quality control and preprocessing of RNA-Seq raw reads
- Normalization of counts
- Differential expression analysis
- Gene ontology and pathway enrichment analysis
- Visualization of results

## Requirements

To run the project, you need to have **R** and the following packages installed:

- **ggplot2** 
- **recount3**
- **recount**
- **edgeR**
- **kableExtra**
- **reshape2**

You can install the required packages with the following command:

```R
install.packages(c("dplyr", "ggplot2", "kableExtra", "reshape2"))
```
and 

```R
install.packages("BiocManager")
BiocManager::install(c("edgeR", "recount", "recount3"))
```





.

