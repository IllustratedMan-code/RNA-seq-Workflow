rdata <- list.files(pattern="\\.Rdata")
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0){
    if(args == "rerun"){
        print("rerunning data")
        r
    }
}
# general library imports ----
library(tximport)
library(readr)
library("pheatmap")

# change these values to your own situations
directory <- "/home/david/Documents/BenoitLab/R/largeRNASEQ"
quantdatafoldername <- "CorrectedQuantData"
fishtype <- "sailfish"

source("/home/david/Documents/BenoitLab/RNA-seq/scripts/functions.r")
txi <- txi(directory, quantdatafoldername, "metadata.csv", fishtype)
Deseq <- initialDeseq(txi[[1]], txi[[2]])

dds <- Deseq[[1]]
res <- Deseq[[2]]
#* vsd ntd rld
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)

ddsMF <- dds
design(ddsMF) <- formula(~part + pesticides)
ddsMF <- DESeq(ddsMF)
Mres <- results(ddsMF)
design(dds) <- formula(~part)
part <- DESeq(dds)
res <- results(part)

DR <- res[ which(Mres$padj < 0.01), ]
Body <- DR[which(DR$log2FoldChange < -2),]
Leg <- DR[which(DR$log2FoldChange > 2), ]

write.table(Body, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/Body.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(Leg, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/Leg.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(res, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/total.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")

directory <- "/home/david/Documents/BenoitLab/RNA-seq"
quantdatafoldername <- "TimeSeriesQuant"
fishtype <- "sailfish"
metadata <- "timeseriesMetadata.csv"
source("/home/david/Documents/BenoitLab/RNA-seq/scripts/functions.r")
txi <- txi(directory, quantdatafoldername, metadata, fishtype)
Deseq <- initialDeseq(txi[[1]], txi[[2]], design = ~Time)
dds <- Deseq[[1]]
res <- Deseq[[2]]
res6 <- results(dds, contrast = c("Time", "6h", "control"))
res2 <- results(dds, contrast = c("Time", "2h", "control"))
res24 <- results(dds, contrast = c("Time", "24h", "control"))
DR6 <- res6[which(res$padj < 0.01), ]
DR2 <- res2[which(res$padj < 0.01), ]
DR24 <- res24[which(res$padj < 0.01), ]

UP6 <- DR6[which(DR6$log2FoldChange > 2), ]
UP2 <- DR2[which(DR2$log2FoldChange > 2), ]
UP24 <- DR24[which(DR24$log2FoldChange > 2), ]

DOWN6 <- DR6[which(DR6$log2FoldChange < -2), ]
DOWN2 <- DR2[which(DR2$log2FoldChange < -2), ]
DOWN24 <- DR24[which(DR24$log2FoldChange < -2), ]

write.table(UP6, file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/UP6.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(UP2, file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/UP2.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(UP24, file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/UP24.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(DOWN6, file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/DOWN6.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(DOWN2, file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/DOWN2.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(DOWN24, file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/DOWN24.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(row.names(UP6), file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/UP6.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(row.names(UP2), file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/UP2.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(row.names(UP24), file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/UP24.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(row.names(DOWN6), file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/DOWN6.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(row.names(DOWN2), file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/DOWN2.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(row.names(DOWN24), file = "/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/timeSeriesGenes/DOWN24.csv", quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")

pestdd <-dds
design(pestdd) <- formula(~part + pesticides)
pestdd <- DESeq(pestdd)
permres <- results(pestdd, contrast = c("pesticides", "perm", "control"))
deetres <- results(pestdd, contrast = c("pesticides", "deet", "control"))

permres <- permres[which(permres$padj < 0.01),]
deetres <- deetres[which(deetres$padj < 0.01),]
upPermvsControl <- permres[which(permres$log2FoldChange >2),]
downPermvsControl <- permres[which(permres$log2FoldChange < -2),]
upDeetvsControl <- deetres[which(deetres$log2FoldChange > 2), ]
downDeetvsControl <- deetres[which(deetres$log2FoldChange < -2), ]
write.table(upPermvsControl, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/upPerm.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(downPermvsControl, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/downPerm.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(upDeetvsControl, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/upDeet.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(downDeetvsControl, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/downDeet.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")

commonUpPermBody = intersect(row.names(Body), row.names(upPermvsControl))
upPermBody = upPermvsControl[commonUpPermBody,]
commonUpPermLeg = intersect(row.names(Leg), row.names(upPermvsControl))
upPermLeg = upPermvsControl[commonUpPermLeg, ]
commonDownPermBody = intersect(row.names(Body), row.names(downPermvsControl))
downPermBody = downPermvsControl[commonDownPermBody,]
commonDownPermLeg = intersect(row.names(Leg), row.names(downPermvsControl))
downPermLeg = downPermvsControl[commonDownPermLeg,]

commonUpDeetBody = intersect(row.names(Body), row.names(upDeetvsControl))
upDeetBody = upDeetvsControl[commonUpDeetBody,]
commonUpDeetLeg = intersect(row.names(Leg), row.names(upDeetvsControl))
upDeetLeg = upDeetvsControl[commonUpDeetLeg,]
commonDownDeetBody = intersect(row.names(Body), row.names(downDeetvsControl))
downDeetBody = downDeetvsControl[commonDownDeetBody,]
commonDownDeetLeg = intersect(row.names(Leg), row.names(downDeetvsControl))
downDeetLeg = downDeetvsControl[commonDownDeetLeg,]


write.table(upPermBody, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/upPermBody.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(upPermLeg, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/upPermLeg.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(upDeetBody, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/upDeetBody.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(upDeetLeg, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/upDeetLeg.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(downPermBody, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/downPermBody.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(downPermLeg, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/downPermLeg.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(downDeetBody, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/downDeetBody.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(downDeetLeg, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/downDeetLeg.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")

newdd <- dds
design(newdd) <- formula(~part + part:pesticides)
newdd <- DESeq(newdd)
newres <- results(newdd)

select <- order(rowMeans(counts(dds, normalized=FALSE)), decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("pesticides","part")])

jpeg("generalfigures/plotMA_deseq.jpg")
plotMA(results)
dev.off()

jpeg("generalfigures/plotcounts_deseq.jpg")
plotCounts(data, gene=which.min(results$padj), intgroup = "part")
dev.off()

jpeg("generalfigures/meanSdvsd_deseq.jpg")
vsn::meanSdPlot(assay(vsd))
dev.off()

jpeg("generalfigures/meadnSdrld_deseq.jpg")
vsn::meanSdPlot(assay(rld))
dev.off()

jpeg("generalfigures/meanSdntd_deseq.jpg")
vsn::meanSdPlot(assay(ntd))
dev.off()

jpeg("generalfigures/pcaplot_deseq.jpg", width = 800, height=800)
plotPCA(vsd, intgroup=c("part", "pesticides"))
dev.off()

jpeg("generalfigures/pheatmapplotntd.jpg", width=800, height=800)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
dev.off()

jpeg("generalfigures/pheatmapplotvsd.jpg", width=800, height=800)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

jpeg("generalfigures/pheatmapplotrld.jpg", width=800, height=800)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

library(tidyverse)
vstarray <- assay(vsd)

range <- vstarray[]
names <- rownames(range)
maxs <- rowMaxs(range)
mins <- rowMins(range)

maxs <- maxs-mins
pt <- tibble(names, maxs)
pt
pt <- pt %>% filter(maxs > 10)
pt
ggplot(pt, aes(x=names, y= maxs)) + geom_col()
descending <- arrange(pt, desc(maxs))
descending <- factor(descending["maxs"])
descending <- tibble(row = seq(length(descending)), descending)
ggplot(descending, aes(x=row, y= descending)) + geom_col()

sum(results$padj < 0.05, na.rm=T)
res <- results

Upregulated <- res[ which(res$padj < 0.03 & res$log2FoldChange > 0), ]
DownRegulated <- res[ which(res$padj < 0.03 & res$log2FoldChange < 0), ]
DownRegulated[order(DownRegulated$padj), ]
head(DownRegulated$padj)
D <- rownames(DownRegulated)
Up <- rownames(Upregulated)
DR <- res[ which(res$padj < 0.03), ]
DR[order(DR$padj),]

write.table(D, file="Genes/DownRegulated.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(Up, file="Genes/UpRegulated.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")

d <- mcols(Deseq,use.names=TRUE)
a <- data.frame(rownumb = seq(nrow(d)), d)
a <- data.frame(rownumb = seq(nrow(d)), basemean = d[,1])
library(ggplot2)
ggplot(a, aes(x=rownumb, y=dispGeneEst)) + geom_point()

Directory <- ("/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/")
DeseqFiles <- list.files(Directory, pattern = "*.csv")
for (i in DeseqFiles) {
  d <- read.csv(paste(Directory, i, sep = ""))
  write.csv(row.names(d), paste(Directory, "genelists/", i, sep = ""), row.names = FALSE, quote = FALSE)
}
