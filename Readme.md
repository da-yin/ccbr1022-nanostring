Nanostring RNA data analysis
=======================

This project is an analysis of Nanostring data on 8 patients who experienced low blood counts after CAR T-cell therapy.



Data source
-----------

The raw data for this project comes from James Kochenderfer Lab at NIH



Analysis
--------

The goals of the analysis are:

- examine Nanostring data from between 2 to 7 time-points for each patient.
- look for differences in gene expression between different time-points.


Analysis code is written for R



---

Introduction
--------

NanoString technology provides a platform to directly measure mRNA expression without cDNA synthesis and amplification steps.The expression level of a target molecule is measured by counting the number of times the barcode for that molecule is detected by a digital analyzer. The system is sensitive enough to detect low abundance molecules. The NanoString nCounter system also provides more accurate quantifications of mRNA expressions than polymerase chain reaction (PCR)-based methods and microarrays in formalin fixed paraffin embedded samples, where RNA degradation is commonly observed. 

This project has a datset with panel of 770 endogenous genes, 10 housekeeping genes, several positive and negative genes. 



Basically, two classes of methods are available to normalize gene expression data using global normalization factors. They are the control-based normalization and the average-bulk normalization. The former class of methods assumes the total expression level summed over a pre-specified group of genes is approximately the same across all the samples. 

The latter class of methods assumes most genes are not significantly Differentially Expressed (DE) across all the samples. The control-based normalization often uses RNA from a group of internal control genes (e.g., housekeeping genes) or external spike-in RNA [e.g., ERCC RNA (Jiang et al., 2011)], while the average-bulk normalization is more commonly used for their universality. 


‘accuracy and reliability of gene expression results are dependent upon the proper normalization of the data against internal reference genes’ (5). Most current normalization methods use spiked-in negative and positive controls and reference (housekeeping) genes, each of which has a good record of effectiveness in related technologies. 

![housekeeping genes](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/housekeepingGenes.pdf)

![normalized gene counts](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/normalizedCounts.jpeg)



