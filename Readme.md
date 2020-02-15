Nanostring RNA data analysis
=======================

This project is an analysis of Nanostring mRNA data on 8 patients who experienced low blood counts after CAR T-cell therapy.



Data source
-----------

The raw data for this project comes from James Kochenderfer Lab at NIH



Analysis
--------

The goals of the analysis are:

- test different normalization methods for scaling the nanostring count data.
- utilize a generalized linear model(GLM) of the negative binomial family to characterize count data.
- examine Nanostring data from between 2 to 7 time-points for each patient.
- look for differences in gene expression between different time-points.


Analysis code is written for R



---

Introduction
--------

NanoString technology provides a platform to directly measure mRNA expression without cDNA synthesis and amplification steps.The expression level of a target molecule is measured by counting the number of times the barcode for that molecule is detected by a digital analyzer. The system is sensitive enough to detect low abundance molecules. The NanoString nCounter system also provides more accurate quantifications of mRNA expressions than polymerase chain reaction (PCR)-based methods and microarrays in formalin fixed paraffin embedded samples, where RNA degradation is commonly observed. 

This project has a datset with panel of 770 endogenous genes, 10 housekeeping genes, several positive and negative genes. For count data, two classes of methods are generally used for normalization. The first one is average-bulk normalization which assumes most genes are not significantly Differentially Expressed (DE) across all the samples. This is commonly used in RNA-seq data normalization such as TMM, UQ, RLE, etc. 

The other class of methods, here used in nanostring data, assumes the total expression level summed over a pre-specified group of genes is approximately the same across all the samples.The control-based normalization often uses RNA from a group of internal control genes (e.g., housekeeping genes) and external spike-in RNA (positive and negative genes). The accuracy and reliability of gene expression results are dependent upon the proper normalization of the data against these internal reference genes.

First I plotted the normalized counts for housekeeping genes. These are ten genes provided by Nanostring, and is expected to have stable and similar level of expression across all samples. However it looks like only G6PD, OAZ1 and SDHA have stable expression. Several other genes even have zeros in some of the samples, which is not expected. 

![housekeeping genes](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/housekeepingGenes.png)

![normalized gene counts](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/normalizedCounts.png)


![normalized gene counts PCA](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/normalizedCounts_PCA.png)


