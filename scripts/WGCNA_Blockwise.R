setwd("/home/david/Documents/BenoitLab/R/largeRNASEQ")
# Packages to install
# install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") )
# source("http://bioconductor.org/biocLite.R")
# biocLite("impute")
# install.packages("WGCNA")
####
# BiocManager::install(c("GO.db","preprocessCore","impute"))
library(WGCNA)
library(flashClust)
# library(tidyr)
options(stringsAsFactors = FALSE)
# enableWGCNAThreads()
allowWGCNAThreads()
# ^Allows WGCNA to run multiple processes at once if you want
exprData <- read.csv("WGCNACorrectedQuant.csv", row.names = 1, header = TRUE)
# gene.names=rownames(exprData)
# ^Input data should be nothing but expression values labelled by sample and gene id
trans.exprData <- t(exprData)
# ^WGCNA expects gene ids to be in columns, this line transposes the data

# Removing genes with 0 variance
gsg <- goodSamplesGenes(trans.exprData, verbose = 3)
gsg$allOK
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing genes:", paste(names(trans.exprData)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing samples:", paste(rownames(trans.exprData)[!gsg$goodSamples], collapse = ", ")))
  }
  # Remove the offending genes and samples from the data:
  trans.exprData <- trans.exprData[gsg$goodSamples, gsg$goodGenes]
}
gene.names <- colnames(trans.exprData)

# Loading in trait data
datTraits <- read.csv("WGCNATruthTable.csv")
# form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) <- datTraits$sample
datTraits$sample <- NULL
table(rownames(datTraits) == rownames(trans.exprData))
# ^should return TRUE if datasets align correctly, otherwise your names are out of order


powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(trans.exprData, dataIsExpr = TRUE, powerVector = powers, corFnc = cor, corOptions = list(use = "p"), networkType = "signed")
setwd("/home/david/Documents/BenoitLab/RNA-seq/figures/WGCNA-Blockwise/")
# Plot results
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
pdf("Connectivity.pdf", width = 5, height = 5)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
# Red line corresponds to using an R^2 cut-off
abline(h = 0.80, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# Generating adjacency and TOM similarity matrices based on the selected softpower
# Choose softpower based on the lowest number that reaches the red line from the Scale independence plot
softPower <- 8
# calclute the adjacency matrix
# adj= adjacency(trans.exprData,type = "signed", power = softPower);

bwnet <- blockwiseModules(trans.exprData,
  maxBlockSize = 16000,
  power = softPower, TOMType = "signed", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = .8,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  saveTOMFileBase = "TOM-blockwise",
  verbose = 3
)

moduleLabels <- bwnet$colors
moduleColors <- labels2colors(bwnet$colors)
MEs <- bwnet$MEs
geneTree <- bwnet$dendrograms[[1]]

pdf("Genetree.pdf", width = 15, height = 5)
plotDendroAndColors(bwnet$dendrograms[[1]], moduleColors[bwnet$blockGenes[[1]]],
  "Module colors",
  dendroLabels = F, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
plotDendroAndColors(bwnet$dendrograms[[2]], moduleColors[bwnet$blockGenes[[2]]],
  "Module colors", #
  dendroLabels = F, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()



# Correlate traits --------------------------------------------------------
# Define number of genes and samples
nGenes <- ncol(trans.exprData)
nSamples <- nrow(trans.exprData)
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(trans.exprData, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


# Print correlation heatmap between modules and traits
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")",
  sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))


# display the corelation values with a heatmap plot
# INCLUE THE NEXT LINE TO SAVE TO FILE
# pdf(file="heatmap.pdf")

pdf("Module relationships.pdf", width = 25, height = 15)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)
dev.off()

write.table(moduleTraitPvalue, "moduleTraitPvalue.txt", sep = "\t")
# For further analysis, if you wanted to pull out genes belonging to a certain module, you can use the following command:
# names(trans.exprData)[dynamicColors=="brown”]


# Creates text files for each color module and a folder to place them in
dir.create("Modules")
module_colors <- setdiff(unique(moduleColors), "")
for (color in module_colors) {
  module <- gene.names[which(moduleColors == color)]
  write.table(module, paste("/home/david/Documents/BenoitLab/RNA-seq/figures/WGCNA-Blockwise/Modules/module_", color, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Heatmap
pdf("heatmap.pdf", width = 25, height = 15)
module.order <- unlist(tapply(1:ncol(trans.exprData), as.factor(moduleColors), I))
m <- t(t(trans.exprData[, module.order]) / apply(trans.exprData[, module.order], 2, max))
heatmap(t(m), zlim = c(0, 1), col = gray.colors(100), Rowv = NA, Colv = NA, labRow = NA, scale = "none", RowSideColors = moduleColors[module.order])
dev.off()

# Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
png(filename = "TOM.png", width = 6, height = 6, units = "in", res = 600)
TOMplot(diss1, hier1, as.character(moduleColors[restGenes]))
dev.off()



# Multi-dimensional scaling plot
pdf("mds plot.pdf", width = 12, height = 12)
cmd1 <- cmdscale(as.dist(dissTOM), 2)
par(mfrow = c(1, 1))
plot(cmd1,
  col = as.character(moduleColors), main = "MDS plot",
  xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2"
)
dev.off()