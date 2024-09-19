library(recount3)
library(recount)
library(edgeR)
library(kableExtra)
library(ggplot2)
library(reshape2)
# import all the data I need: one referring to the brain, one referring to the lung and one referring to the pancreas
# these data are chosen "randomly", referring to the fifth and sixth digits of the university number: the brain is by default and the others are 5.6
rse_brain <- readRDS("rse_brain.RDS") 
rse_lung <- readRDS("rse_lung.RDS")
rse_pancreas <- readRDS("rse_pancreas.RDS")

assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_lung)$counts <- transform_counts(rse_lung)
assays(rse_pancreas)$counts <- transform_counts(rse_pancreas)

# select 3 samples for each tissues, also here the samples' number is randomized using the first and the second number of the 
# university ID number, in my case 45. If the sample selected don't pass these "test":

# Observing the minimun RIN
# colData(rse_pancreas)$gtex.smrin[45]
# colData(rse_brain)$gtex.smrin[46]         # if the output is > 6 I can use the sample, otherwise I add +1 and choose
# colData(rse_brain)$gtex.smrin[47]         # another sample

# Estimated fraction of rRNA
# colData(rse_pancreas)$gtex.smrrnart[45]
# colData(rse_brain)$gtex.smrrnart[46]      # if the output is < 10% I can use the sample, otherwise I add +1 and choose
# colData(rse_brain)$gtex.smrrnart[47]      # another sample

# The percentage of mapped reads
# colData(rse_pancreas)$"recount_qc.star.uniquely_mapped_reads_%_both"[48]
# colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[46]     # if the output is > 85% I can use the sample,
# colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[47]     # otherwise i add +1 and choose another sample


x11()
par(mfrow = c(3,1))
rinplot <- colData(rse_brain)$gtex.smrin
boxplot(rinplot,
        horizontal = TRUE,
        main = "RIN")
boxplot(colData(rse_brain)$gtex.smrrnart,
        horizontal = TRUE,
        main = "rRNA fraction")
mapped1 <- colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped1, 
        horizontal = TRUE,
        main = "% of uniquely mapping")

x11()
par(mfrow = c(3,1))
rinplot2 <- colData(rse_lung)$gtex.smrin
boxplot(rinplot2,
        horizontal = TRUE,
        main = "RIN")
boxplot(colData(rse_lung)$gtex.smrrnart,
        horizontal = TRUE,
        main = "rRNA fraction")
mapped2 <- colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped2, 
        horizontal = TRUE,
        main = "% of uniquely mapping")

x11()
par(mfrow = c(3,1))
rinplot3 <- colData(rse_pancreas)$gtex.smrin
boxplot(rinplot3,
        horizontal = TRUE,
        main = "RIN")
boxplot(colData(rse_pancreas)$gtex.smrrnart,
        horizontal = TRUE,
        main = "rRNA fraction")
mapped3 <- colData(rse_pancreas)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped3, 
        horizontal = TRUE,
        main = "% of uniquely mapping")

rse_brain_selected <- rse_brain[,c(48,49,50)]
rse_lung_selected <- rse_lung[,c(49,46,47)]
rse_pancreas_selected <- rse_pancreas[,c(45,48,47)]
counts_brain_selected <- assays(rse_brain_selected)$counts
counts_lung_selected <- assays(rse_lung_selected)$counts
counts_pancreas_selected <- assays(rse_pancreas_selected)$counts

# create a data frame with the selected data and also arrange row names and column names
x <- cbind(counts_brain_selected,counts_lung_selected,counts_pancreas_selected)

# as column put the samples name
colnames(x) <- c("Brain48", "Brain49","Brain50","Lung49", "Lung46","Lung47","Pancreas45","Pancreas48","Pancreas47")

# as row put genes name
rownames(x) <- rowData(rse_brain_selected)$gene_name

# DGEList: This function creates an object of class DGEList. "DGE" stands for "Differential Gene Expression". 
# This object is a specialized data structure that edgeR uses to store and manipulate gene count data.
y <- DGEList(counts=x)

# defie the belonging group for each samples 
group <- as.factor(c("Brain","Brain","Brain","Lung","Lung","Lung","Pancreas","Pancreas","Pancreas"))

y$samples$group <- group

# adding some information as column which can be useful for the study/analysis
y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,
                             colData(rse_lung_selected)$gtex.smrin,
                             colData(rse_pancreas_selected)$gtex.smrin))

y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,
                               colData(rse_lung_selected)$gtex.smtsd,
                               colData(rse_pancreas_selected)$gtex.smtsd))

y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,
                             colData(rse_lung_selected)$gtex.sex,
                             colData(rse_pancreas_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,
                             colData(rse_lung_selected)$gtex.age,
                             colData(rse_pancreas_selected)$gtex.age))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,
                              colData(rse_lung_selected)$gtex.smrrnart,
                              colData(rse_pancreas_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", 
                                colData(rse_lung_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",
                                colData(rse_pancreas_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", 
                              colData(rse_lung_selected)$"recount_qc.aligned_reads%.chrm",
                              colData(rse_pancreas_selected)$"recount_qc.aligned_reads%.chrm"))
y

# check how many genes have zero counts, and then we use the edgeR function for removing zero or 
# low expression genes with default parameters (notice that we have to specify to it how samples are split into groups).
table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

dim(y)

# cpm function: The cpm function in edgeR calculates counts per million (CPM) for the gene expression data. 
# CPM is a normalization method that scales the raw gene counts to account for differences in sequencing depth 
# between samples, making the data more comparable across samples.
# log=TRUE Argument: The log=TRUE argument specifies that the counts per million should be log-transformed. 
# By default, a log2 transformation is applied, which means the function calculates log2(CPM + 1) 
# to avoid taking the log of zero.
logcpm_before <- cpm(y, log=TRUE)

# calcNormFactors Function: This function calculates normalization factors to scale the raw library sizes in a way 
# that makes the expression levels more comparable between samples. These factors account for differences in sequencing 
# depth and RNA composition between the samples.
# method = "TMM" Argument: The method argument specifies the normalization method to be used. "TMM" stands for "Trimmed Mean of M-values". 
# TMM is a widely used method for normalizing RNA-Seq data that accounts for differences in library sizes and compositional
# biases. It trims the most extreme log-fold changes and the most lowly expressed genes, then computes a weighted average of the 
# remaining log-fold changes to estimate the relative scaling factors between libraries.
y <- calcNormFactors(y, method = "TMM")
y

# I run the same function, but now the data are normalized, in fact I can see differences in the box-plots.
logcpm_after <- cpm(y, log=TRUE)

boxplot(logcpm_before, notch=T, ylab = "Log CPM Before Normalization", xlab = "Tissue samples")
boxplot(logcpm_after, notch=T, ylab = "Log CPM After Normalization", xlab = "Tissue samples")

normalization_data <- data.frame(Sample = y$samples$group, NormalizationFactor = y$samples$norm.factors)
normalization_kable = kbl(normalization_data, align = "c") %>%
  kable_material()

# model.matrix Function: This function creates a design matrix for linear modeling. The design matrix is used to specify 
# the experimental design and the relationships between the samples and the conditions being tested.
# ~0+group Formula: This formula specifies the structure of the design matrix:

# ~ indicates that the right side of the formula specifies the model.
# 0+ indicates that the model should not include an intercept term. By omitting the intercept (using 0+), 
# the design matrix will have a separate column for each level of the group factor, instead of having an intercept
# and one fewer column.

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

# plotMDS function: This function creates a multidimensional scaling (MDS) plot. An MDS plot is a type of unsupervised 
# analysis that visualizes the similarity or dissimilarity of data points (in this case, RNA-Seq samples) in a reduced 
# dimensional space. It is particularly useful for exploring patterns and clusters in high-dimensional data, like gene 
# expression data.
x11()
par(mfrow=c(1,5))
plotMDS(logcpm_after, labels=group, main = "Histology")
plotMDS(logcpm_after, labels=y$samples$rRNA, main = "rRNA fraction")
plotMDS(logcpm_after, labels=y$samples$chrm, main = "% mitochondrial RNA")
plotMDS(logcpm_after, labels=y$samples$age, mai = "Age")
plotMDS(logcpm_after, labels=y$samples$sex, mai = "Sex")


# estimateDisp Function: This function estimates the dispersion parameters for the RNA-Seq data. Dispersion in RNA-Seq data 
# refers to the variability in gene expression that cannot be explained by the experimental design alone. 
# Accurate estimation of dispersion is crucial for reliable identification of differentially expressed genes.
y <- estimateDisp(y, design)

# plotBCV Function: This function generates a Biological Coefficient of Variation (BCV) plot. The BCV plot is a graphical 
# representation that shows the estimated dispersions (variability) of the gene expression data. It helps in visualizing 
# the relationship between gene expression levels and their variability across samples.
plotBCV(y)

# glmQLFit Function: This function fits a quasi-likelihood (QL) negative binomial generalized linear model (GLM) to the 
# RNA-Seq data. The quasi-likelihood approach accounts for variability in the dispersion estimates, providing more robust 
# statistical inference compared to traditional negative binomial GLM.
fit <- glmQLFit(y, design)

#lung (top) vs brain (bottom)
qlfLB <- glmQLFTest(fit, contrast=c(-1,1,0))
#pancreas (top) vs brain (bottom)
qlfPB <- glmQLFTest(fit, contrast=c(-1,0,1))
#pancreas (top) vs lung (bottom)
qlfPL <- glmQLFTest(fit, contrast=c(0,-1,1))

qlfLB
qlfPB
qlfPL

head(qlfLB$table)
head(qlfPB$table)
head(qlfPL$table)

up_down_counts <- function(qlf, threshold = 0.01) {
  upregulated <- sum(qlf$table$PValue < threshold & qlf$table$logFC > 0)
  downregulated <- sum(qlf$table$PValue < threshold & qlf$table$logFC < 0)
  return(c(Upregulated = upregulated, Downregulated = downregulated))
}

counts_LB <- up_down_counts(qlfLB)
counts_PB <- up_down_counts(qlfPB)
counts_PL <- up_down_counts(qlfPL)

comparison_names <- c("Lung vs Brain", "Pancreas vs Brain", "Pancreas vs Lung")
upregulated_counts <- c(counts_LB["Upregulated"], counts_PB["Upregulated"], counts_PL["Upregulated"])
downregulated_counts <- c(counts_LB["Downregulated"], counts_PB["Downregulated"], counts_PL["Downregulated"])

results_df <- data.frame(
  Comparison = comparison_names,
  Upregulated = upregulated_counts,
  Downregulated = downregulated_counts
)

results_melted <- melt(results_df, id.vars = "Comparison")

ggplot(results_melted, aes(x = Comparison, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("white","grey")) +
  labs(y = "Number of Genes", x = "Comparison", fill = "Regulation") +
  theme_classic()+
  theme(legend.position = "none")

kable_up_down = kbl(results_df, align = "c") %>%
  kable_material()


# topTags Function: This function extracts and summarizes the top differentially expressed genes from the results 
# of a statistical test. It provides a ranked table of genes based on the specified criteria.
# adjust.method = "BH" Argument: The adjust.method argument specifies the method for adjusting p-values to account for 
# multiple testing. "BH" stands for the Benjamini-Hochberg method, which controls the false discovery rate (FDR). 
# This method is commonly used in RNA-Seq data analysis to reduce the likelihood of false positives.
topTags(qlfLB, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfPB, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfPL, n=10,adjust.method = "BH", sort.by = "PValue")

# save the information of topTags of all genes for each tissue
resultsLB <- topTags(qlfLB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsPB <- topTags(qlfPB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsPL <- topTags(qlfPL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

write.table(resultsLB, "resultsLB.txt")
write.table(resultsPB, "resultsPB.txt")
write.table(resultsPL, "resultsPL.txt")

summary(decideTests(qlfLB, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfLB, p.value=0.05, adjust.method = "BH", lfc=1))
summary(decideTests(qlfLB, p.value=0.01, adjust.method = "BH", lfc=0))
summary(decideTests(qlfLB, p.value=0.01, adjust.method = "BH", lfc=1))

# Genes up-regulated lung vs brain
LBtest <- decideTests(qlfLB, p.value=0.01, adjust.method = "BH",lfc=1)
FDR_compliant <- rownames(resultsLB[resultsLB@.Data[[1]]$FDR < 0.01,])

upregulated_LB <- rownames(LBtest)[LBtest == -1]
upregulated_LB_FDR_compliant <- intersect(upregulated_LB, FDR_compliant)

# Genes up-regulated pancreas vs brain
PBtest <- decideTests(qlfPB, p.value=0.01, adjust.method = "BH",lfc=1)
FDR_compliant <- rownames(resultsPB[resultsPB@.Data[[1]]$FDR < 0.01,])

upregulated_PB <- rownames(PBtest)[PBtest == -1]
upregulated_PB_FDR_compliant <- intersect(upregulated_PB, FDR_compliant)

# Genes up-regulated in both
brain_upregulated <- intersect(upregulated_LB_FDR_compliant,upregulated_PB_FDR_compliant)
length(brain_upregulated) # 2588

write.table(brain_upregulated, "Brain_genes.txt")

PBtest <- decideTests(qlfPB, p.value=0.01, adjust.method = "BH",lfc=1)
FDR_compliant <- rownames(resultsPB[resultsPB@.Data[[1]]$FDR < 0.01,])

upregulated_PB2 <- rownames(PBtest)[PBtest == +1]
upregulated_PB_FDR_compliant2 <- intersect(upregulated_PB2, FDR_compliant)

PLtest <- decideTests(qlfPL, p.value=0.01, adjust.method = "BH",lfc=1)
FDR_compliant <- rownames(resultsPL[resultsPL@.Data[[1]]$FDR < 0.01,])

upregulated_PL <- rownames(PLtest)[PLtest == +1]
upregulated_PL_FDR_compliant <- intersect(upregulated_PL, FDR_compliant)

# Genes up-regulated in both
pancreas_upregulated <- intersect(upregulated_PL_FDR_compliant,upregulated_PB_FDR_compliant2)
length(pancreas_upregulated) # 906

write.table(pancreas_upregulated, "Pancreas_genes.txt")

LBtest <- decideTests(qlfLB, p.value=0.01, adjust.method = "BH",lfc=1)
FDR_compliant <- rownames(resultsLB[resultsLB@.Data[[1]]$FDR < 0.01,])

upregulated_LB2 <- rownames(LBtest)[LBtest == +1]
upregulated_LB_FDR_compliant2 <- intersect(upregulated_LB2, FDR_compliant)

PLtest <- decideTests(qlfPL, p.value=0.01, adjust.method = "BH",lfc=1)
FDR_compliant <- rownames(resultsPL[resultsPL@.Data[[1]]$FDR < 0.01,])

upregulated_PL2 <- rownames(PLtest)[PLtest == -1]
upregulated_PL_FDR_compliant2 <- intersect(upregulated_PL2, FDR_compliant)

# Genes up-regulated in both
lung_upregulated <- intersect(upregulated_PL_FDR_compliant2,upregulated_LB_FDR_compliant2)
length(lung_upregulated) # 1589

write.table(lung_upregulated, "Lung_genes.txt")

up_dataframe = data.frame(
  Tissue = c("Brain", "Lung", "Pancreas"),
  Upregulated = c(length(brain_upregulated), length(lung_upregulated), length(pancreas_upregulated))
)

ggplot(up_dataframe, aes(x = Tissue, y = Upregulated, fill = Tissue)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("darkgrey","lightgrey","white")) +
  labs(y = "Number of Uniquely Upregulated Genes", x = "Tissue") +
  theme_classic() +
  theme(legend.position = "none")

up_kable = kbl(up_dataframe, align = "c") %>%
  kable_material()

assays(rse_brain)$TPM <- getTPM(rse_brain)
assays(rse_lung)$TPM <- getTPM(rse_lung)
assays(rse_pancreas)$TPM <- getTPM(rse_pancreas)
which(rowData(rse_brain)$gene_name == "INS")

boxplot(assays(rse_brain)$TPM[12801,],assays(rse_lung)$TPM[12801,], assays(rse_pancreas)$TPM[12801,], outline=F,
        names = c("Brain", "Lung", "Pancreas"),
        ylab = "Expression in TPM",
        main = "Gene: INS")

which(rowData(rse_brain)$gene_name == "GRIN1")

boxplot(assays(rse_brain)$TPM[51767,],assays(rse_lung)$TPM[51767,], assays(rse_pancreas)$TPM[51767,], outline=F,
        names = c("Brain", "Lung", "Pancreas"),
        ylab = "Expression in TPM",
        main = "Gene: GRIN1")

which(rowData(rse_brain)$gene_name == "SFTPC")

boxplot(assays(rse_brain)$TPM[48273,],assays(rse_lung)$TPM[48273,], assays(rse_pancreas)$TPM[48273,], outline=F,
        names = c("Brain", "Lung", "Pancreas"),
        ylab = "Expression in TPM",
        main = "Gene: SFTPC")

# test upregulated brain GRIN1
wilcox.test(assays(rse_brain)$TPM[51767,],assays(rse_lung)$TPM[51767,], paired = F, alternative = "greater")
wilcox.test(assays(rse_brain)$TPM[51767,],assays(rse_pancreas)$TPM[51767,], paired = F, alternative = "greater")
# p.value < 2.2e-16 reject Ho 

# test upregulated lung SFTPC
wilcox.test(assays(rse_lung)$TPM[48273,],assays(rse_brain)$TPM[48273,], paired = F, alternative = "greater")
wilcox.test(assays(rse_lung)$TPM[48273,],assays(rse_pancreas)$TPM[48273,], paired = F, alternative = "greater")
# p.value < 2.2e-16 reject Ho 

wilcox.test(assays(rse_pancreas)$TPM[12801,],assays(rse_brain)$TPM[12801,], paired = F, alternative = "greater")
wilcox.test(assays(rse_pancreas)$TPM[12801,],assays(rse_lung)$TPM[12801,], paired = F, alternative = "greater")
# p.value < 2.2e-16 reject Ho 
