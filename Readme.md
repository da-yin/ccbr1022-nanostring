Nanostring RNA data analysis
=======================

This project is an analysis of Nanostring mRNA data on 8 patients who experienced low blood counts after CAR T-cell therapy.


Data source
-----------

The raw data for this project comes from James Kochenderfer Lab at NIH
- /RawData/Cytopenia_nanostring has the nanostring RCC files, each represent a sample
- /RawData/ccbr1022_rawcounts.csv has the raw counts parsed from RCC combining all samples.
- /RawData/ccbr1022_metadata.csv has the metadata for all samples including treatment, cytopenia severity, patient number, sample number etc. 


Objectives
--------

- test different normalization methods for scaling the nanostring count data.
- utilize a generalized linear model(GLM) of the negative binomial family to characterize count data.
- examine Nanostring data from between 2 to 7 time-points for each patient.
- look for differences in gene expression between different time-points.


Analysis code is written for R version 3.6.1

---

Summary
--------

NanoString technology provides a platform to directly measure mRNA expression without cDNA synthesis and amplification steps.The expression level of a target molecule is measured by counting the number of times the barcode for that molecule is detected by a digital analyzer. The system is sensitive enough to detect low abundance molecules. The NanoString nCounter system also provides more accurate quantifications of mRNA expressions than polymerase chain reaction (PCR)-based methods and microarrays in formalin fixed paraffin embedded samples, where RNA degradation is commonly observed. 

This project has a datset with panel of 770 endogenous genes, 10 housekeeping genes, several positive and negative genes. For count data, two classes of methods are generally used for normalization. The first one is average-bulk normalization which assumes most genes are not significantly Differentially Expressed (DE) across all the samples. This is commonly used in RNA-seq data normalization such as TMM, UQ, RLE, etc. 

The other class of methods, here used in nanostring data, assumes the total expression level summed over a pre-specified group of genes is approximately the same across all the samples.The control-based normalization often uses RNA from a group of internal control genes (e.g., housekeeping genes) and external spike-in RNA (positive and negative genes). The accuracy and reliability of gene expression results are dependent upon the proper normalization of the data against these internal reference genes.

First I plotted the normalized counts for housekeeping genes **Figure1**. These are ten genes provided by Nanostring, and is expected to have stable and similar level of expression across all samples. However it looks like only G6PD, OAZ1 and SDHA have stable expression. Several other genes even have zeros in some of the samples, which is not expected. 

Then I used use the reference genes and normalized the counts **Figure2**. Then I used PCA to show the principal causes of variation in a dataset **Figure3**. Then 200 most variable genes are used to draw heatmap and show unsupervised clustering of samples **Figure4**. MD plot using Nanostring provided housekeeping genes shows that many of the differentilly expressed genes are downregulated (>80%) **Figure5**. 

I used DESeq2 to calculate differentilly expressed genes between the cytopenia and none-cytopenia samples. Then I selected the genes with the 100th largest P values (most constant genes). Then out of these genes, the ones that have an above median average total counts are selected as new reference housekeeping genes for subsequent NanostringDIff analysis. **Figure6** shows the fold change vs. mean expression. **Figure7** shows the DEGs calcualted using nanostringDiff and new reference genes and clustering heatmap.

### Figure1 housekeeping genes
![housekeeping genes](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/housekeeping_normalized.png)

### Figure2 normalized gene counts
![normalized gene counts](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/normalizedCounts.png)

### Figure3 normalized gene counts PCA
![normalized gene counts PCA](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/PCA_normalized_batchCorr.png)

### Figure4 unsupervised clustering of samples based on 200 most variable genes
![200 most variable genes heatmap](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/Heatmap_normalized_batchcor.png)

### Figure5 MD plot using nanostring provided housekeeping for differential expression
![MD before](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/MDplot_before.png)

### Figure6 MD plot using DESeq2 selected housekeeping for differential expression
![MD after](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/MDplot_after.png)

### Figure7 supervised clustering of samples based on DEGs
![DEGs genes heatmap](https://github.com/da-yin/ccbr1022-nanostring/blob/master/Analysis/Results/DEGheatmap_supervised.png)





