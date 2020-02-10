suppressMessages(library(Biobase))
suppressMessages(library(NanoStringDiff))
suppressMessages(library(enrichplot))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(clusterProfiler))
suppressMessages(library(argparse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(ReactomePA))
suppressMessages(library(ggplot2))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ggfortify))
suppressMessages(library(ggrepel))
suppressMessages(library(VennDiagram))
suppressMessages(library(biomaRt))
suppressMessages(library(enrichplot))
suppressMessages(library(gridExtra))
suppressMessages(library(gplots))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(Glimma))
suppressMessages(library(ggpubr))

###################################NanoStringDiff######################################################

setwd("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/")
rawcounts = read.csv(file = "~/Desktop/active_projects/ccbr1022-nanostring/RawData/ccbr1022_rawcounts.csv", header = T)
colnames(rawcounts) = gsub("X","",colnames(rawcounts))

severity = as.factor(str_sub(colnames(rawcounts)[4:ncol(rawcounts)], -4,-1))
group = as.factor(ifelse(severity=="none", "notCytopenic", "cytopenic"))
patient = as.factor(str_sub(colnames(rawcounts)[4:ncol(rawcounts)], 6,7))

designs = data.frame(group, patient)

NanoStringData <- createNanoStringSetFromCsv("~/Desktop/active_projects/ccbr1022-nanostring/RawData/ccbr1022_rawcounts.csv",
                                             header = TRUE, designs)

#normalize raw data using positive, negative, and housekeeping factors
NanoStringData=estNormalizationFactors(NanoStringData)

c = positiveFactor(NanoStringData)
d = housekeepingFactor(NanoStringData)
k = c * d
lamda_i = negativeFactor(NanoStringData)
Y = exprs(NanoStringData)
Y_n = sweep(Y, 2, lamda_i, FUN = "-")
Y_nph = sweep(Y_n, 2, k, FUN = "/")
colnames(Y_nph) = gsub("X","",colnames(Y_nph))
Y_nph[Y_nph <= 0] = 1
Y_nph_log <- log(Y_nph)
Y_nph_palantir = Y_nph

#################################examine the normalized raw data##################################################

par(mar=c(10.1,4.1,4.1,2.1))
long_Y_nph = melt(Y_nph, value.name = "normalized_Counts", varnames=c('gene', 'sample'))
long_Y_nph$cytopenia = str_sub(long_Y_nph$sample, -4,-1)
long_Y_nph$cytopenia=as.factor(long_Y_nph$cytopenia)

#png(filename = "~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/normalizedCounts.png")
pdf("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/normalizedCounts.pdf")
ggplot(long_Y_nph, aes(x=sample, y=normalized_Counts, fill=cytopenia)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(0, 300)
dev.off()

#################################examine sample clusterings##################################################

pca_plot = function(data){
  counts = data.frame(t(data))
  counts$group = str_sub(rownames(counts), -4,-1)
  df = counts[,1:ncol(counts)-1]
  df = df[,apply(df,2,var)!= 0]
  autoplot(prcomp(df, scale. = T), data = counts, colour= 'group',label=T,repel = TRUE,
           label.size = 3,label.repel=T)
}

pdf("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/normalizedCounts_PCA.pdf")
pca_plot(Y_nph)
dev.off()

#scaling and normalizing the data to get counts per thousand
for (c in colnames(counts)){
  if (is.integer(counts[[c]])){
    counts[[c]] = (counts[[c]] * 1000)/sum(counts[[c]])
  }
}

#################################examine housekeeping genes##################################################

#checking to see if the housekeeping genes actually have consistant expression across samples
counts_housekeeping = counts[which(counts$Code.Class=="Housekeeping"),]
counts_housekeeping_tidy = melt(data = counts_housekeeping, value.name = "CPK", variable.name = "sample")
counts_housekeeping_tidy$cytopenia = str_sub(counts_housekeeping_tidy$sample, -4,-1)
housekeeping.genes = unique(counts_housekeeping_tidy$Name)

PlotHouse = function(gene){
  df = counts_housekeeping_tidy[which(counts_housekeeping_tidy$Name==deparse(substitute(gene))),]
  df_CPK = df$CPK
  ggplot(data = counts_housekeeping_tidy[which(counts_housekeeping_tidy$Name==deparse(substitute(gene))),], 
         aes(x = sample, y=CPK, fill=cytopenia)) + geom_bar(stat="identity")+
    theme(axis.text.x = element_text(angle = 0, hjust = 1))+ labs(title = deparse(substitute(gene)))+
    geom_hline(yintercept = mean(df_CPK), color="blue")+geom_hline(yintercept = sd(df_CPK)+mean(df_CPK), color="red")+
    geom_hline(yintercept =mean(df_CPK)-sd(df_CPK), color="red")
}

pdf("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/housekeepingGenes.pdf")
ggarrange(PlotHouse(ABCF1),PlotHouse(G6PD),PlotHouse(NRDE2),PlotHouse(OAZ1),PlotHouse(POLR2A),
          PlotHouse(SDHA),PlotHouse(STK11IP),PlotHouse(TBC1D10B),PlotHouse(TBP),PlotHouse(UBB),
          nrow = 5, ncol = 2)
dev.off()


##############################use DESeq2 to choose new housekeeping genes#############################################
require(DESeq2)

DESeq2rawcounts = read.csv("~/Desktop/active_projects/ccbr1022-nanostring/RawData/ccbr1022_rawcounts.csv", header = T,row.names = 2)
metadata = read.csv("~/Desktop/active_projects/ccbr1022-nanostring/RawData/ccbr1022_metadata.csv", header = T)

DESeq2rawcounts$Code.Class=NULL
DESeq2rawcounts$Accession=NULL
colnames(DESeq2rawcounts)=gsub("X","",colnames(DESeq2rawcounts))
metadata$filename=NULL
#It is critical columns of the count matrix and the rows of the column data (information about samples) are in the same order
rownames(metadata) = metadata$sampleName
colnames(metadata)
all(rownames(metadata) == colnames(DESeq2rawcounts))

dds = DESeqDataSetFromMatrix(countData = DESeq2rawcounts,
                             colData = metadata,
                             design = ~ cytopenia + patient)

# filter out low count genes, total counts for each gene must be equal to lower quantile or higher
keep = rowSums(counts(dds)) >= quantile(rowSums(counts(dds)), probs=c(0.25))
dds = dds[keep,]
dds$condition <- relevel(dds$cytopenia, ref = "NotCytopenic")

dds = DESeq(dds)
res = results(dds)
# choose genes that has the largest P value, and above 50% expression level as new housekeeping genes
constant.genes = rownames(res[order(res$pvalue, decreasing = T),][1:100,])
constant.genes = intersect(constant.genes, rownames(DESeq2rawcounts)[rowSums(DESeq2rawcounts) >= quantile(rowSums(DESeq2rawcounts), probs=c(0.5))])

write.csv(res, "~/Desktop/active_projects/ccbr1022-nanostring/Analysis/processedData/DESeq2.res.csv",
          quote = F, row.names = F)

###swap out the default housekeeping genes with new calcualted housekeeping genes
rawcounts2 = read.csv("~/Desktop/active_projects/ccbr1022-nanostring/RawData/ccbr1022_rawcounts.csv", header = T)
colnames(rawcounts2)=gsub("X","",colnames(rawcounts2))
new.housekeeping = rawcounts2[(rawcounts2$Name%in%constant.genes),]

summary(rawcounts2$Code.Class)
new.housekeeping$Code.Class=as.factor("Housekeeping")
rawcounts2 = rawcounts2[which(rawcounts2$Code.Class!="Housekeeping"),]
rawcounts3= rbind(new.housekeeping,rawcounts2)
summary(rawcounts3$Code.Class)
write.csv(rawcounts3, "~/Desktop/active_projects/ccbr1022-nanostring/Analysis/processedData/ccbr1022_rawcounts_newHousekeeping.csv",
          quote = F, row.names = F)

##############################use nanostringdiff to find DEGs#############################################

rawcounts = read.csv(file = "~/Desktop/active_projects/ccbr1022-nanostring/Analysis/processedData/ccbr1022_rawcounts_newHousekeeping.csv", header = T)
colnames(rawcounts) = gsub("X","",colnames(rawcounts))

severity = as.factor(str_sub(colnames(rawcounts)[4:ncol(rawcounts)], -4,-1))
group = as.factor(ifelse(severity=="none", "notCytopenic", "cytopenic"))
patient = as.factor(str_sub(colnames(rawcounts)[4:ncol(rawcounts)], 6,7))

designs = data.frame(group, patient)

NanoStringData <- createNanoStringSetFromCsv("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/processedData/ccbr1022_rawcounts_newHousekeeping.csv",
                                             header = TRUE, designs)
pheno = pData(NanoStringData)[, c("group","patient")]
group = pheno$group
patient = pheno$patient
design.full = model.matrix(~group+patient)

#normalize raw data using positive, negative, and housekeeping factors
NanoStringData=estNormalizationFactors(NanoStringData)

Result_cytopenia.vs.noCytopenia = glm.LRT(NanoStringData,design.full, Beta = 2)
Result_cytopenia.vs.noCytopenia$table

#Because code automatically set "Cytopenia" as the baseline level of your covariant "group" 
#It will not affect significant test, p-value and q-value. But for logFC, calculate it as -logFC of the output.
Result_cytopenia.vs.noCytopenia$table$logFC=Result_cytopenia.vs.noCytopenia$table$logFC * (-1)

write.csv(Result_cytopenia.vs.noCytopenia$table, "~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/Result_cytopenia.vs.noCytopenia.csv")

###########################generate MD plots###############################################################

DEG = Result_cytopenia.vs.noCytopenia$table
class(DEG)
DEG$sig = ifelse(abs(DEG$logFC)>1 & DEG$qvalue<0.05, "sig","ns" )
DEG$gene = rownames(DEG)

normalized_log = data.frame(Y_nph_log)
normalized_log$ave_log = apply(normalized_log, 1, mean)
normalized_log$gene = rownames(normalized_log)
merged = merge(DEG, normalized_log, by="gene")

pdf("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/MDplot.pdf")
ggplot(merged, aes(x=merged$ave_log, y=merged$logFC, color=merged$sig, ,label=T,repel = TRUE))+ geom_point() + ylim(-40,40)+
  xlab("log mean expression") + ylab("log fold change") + labs(title ="MD plot", col="DEG")+theme(axis.text=element_text(size=12),
                                                                                                  axis.title=element_text(size=14))
dev.off()
################generate heatmap of top DEG genes###############
normalized = data.frame(Y_nph)
colnames(normalized) = gsub("X","",colnames(normalized))
normalized$ave = apply(normalized, 1, mean)
normalized$gene = rownames(normalized)

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

merged2 = merge(DEG, normalized, by="gene")
merged2_sig = merged2[which(merged2$sig=="sig"),]
merged2_sig = merged2_sig[order(merged2_sig$qvalue, decreasing = F),]
rownames(merged2_sig)=merged2_sig$gene
merged2_sig_counts = subset(merged2_sig, select=-c(gene,logFC,lr,pvalue,qvalue,sig,ave))
pdf("~/Desktop/active_projects/ccbr1022-nanostring/Analysis/Results/DEGheatmap.pdf")
heatmap.2(data.matrix(merged2_sig_counts),col=rev(morecols(50)),trace="none", 
          main="DEG genes across samples",scale="row",srtCol=45,cexCol=1,cexRow = 0.7)
dev.off()

##############save Rdata##################
save.image(file = "nanostring")
