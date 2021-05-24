vignette(package="dplyr")

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
design(ddsMF) <- formula(~pesticides)
ddsMF <- DESeq(ddsMF)
Mres <- results(ddsMF)
design(dds) <- formula(~part)
part <- DESeq(dds)
res <- results(part)

DR <- Mres[ which(Mres$padj < 0.01), ]
Body <- DR[which(DR$log2FoldChange < -2),]
Leg <- DR[which(DR$log2FoldChange > 2), ]

write.table(Body, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/Body.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(Leg, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/Leg.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(res, file="/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/total.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")

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

#! EdgeR section ----
library(edgeR)
# EdgR dataset from tximport
counts <- txi$counts
length <- txi$length
normMat <- length / exp(rowMeans(log(length)))
normCts <- counts / normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
y <- DGEList(counts)
y <- scaleOffset(y, normMat)



# EdgR heatmap
logcpm <- cpm(EdgeR)
meta <- read.csv("metadata.csv", sep= ",")
pest <- factor(meta$pesticides, levels=c("perm", "deet", "control"))
part <- factor(meta$part, levels=c("leg", "body"))
design <- model.matrix(~0 + pest + part)
design
keep <- filterByExpr(y, design)
EdgeR <- y[keep, ]
estimateDisp(EdgeR)
y <- estimateGLMCommonDisp(EdgeR, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
est <- estimateDisp(EdgeR, design, robust=TRUE)
design
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit, coef=1:4)
test <- glmQLFTest(fit, coef=1:4)
lrc$PValue
lrt

test$logFC.pestperm - test$logFC.pestdeet

ggplot(test$table, aes(y=logFC.partleg, x=seq(nrow(test)))) + geom_hex()
#additional analysis----
table <- test$table
body <- table[table$logFC.partbody > 2, ]
leg <- table[table$logFC.partbody < -2, ]
table <- body
diff <- table$logFC.pestcontrol - table$logFC.pestdeet
library(tidyverse)
data <- tibble(data = diff, rows = seq(length(diff)))
ggplot(data, aes(x=rows, y=data)) + geom_point()
table$diff <- diff
cont <- table[table$diff > 2, ]
deet <- table[table$diff < -2, ]
diff2 <- table$logFC.pestcontrol - table$logFC.pestperm
table$diff2 <- diff2
cont2 <- table[table$diff2 > 2, ]
perm <- table[table$diff2 < -2, ]



library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
Legcont1 <- cont
Legcont2 <- cont2
deet1 <- deet
perm1 <- perm
d <- data.frame(Pest = c(rep(c("DEET", "DEET", "PERM", "PERM"))),
                Part = c(rep("Body", 4), rep("Leg", 4)),
                Names = c(rep(c("Up Regulated", "Down Regulated"), 4)),
                GeneCount = c(nrow(deet), nrow(cont), nrow(perm), nrow(cont2), nrow(deet1), nrow(Legcont1), nrow(perm1), nrow(Legcont2)))

write.table(d, file="Genes/regulation.csv", quote=FALSE, col.names=TRUE, row.names=FALSE, sep=",")
write.table(rownames(perm), file="Genes/perm.csv", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")

pl ggplot(d, aes(x=Names, y = GeneCount, fill=Names)) + geom_col()
venn.diagram(
  x = list(rownames(deet), rownames(cont)),
  category.names = c("DEET" , "Control"),
  filename = "Deetvcontorlvleg.png",
  compression = "lzw",


  lwd = 2,
  lty = 'blank',
  fill = c("red", "blue"),
)
grid.draw(tmp)

write.table(rownames(perm), file="Genes/perm.csv", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(rownames(deet), file="Genes/deet.csv", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(rownames(cont), file="Genes/deetcont.csv", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(rownames(cont2), file="Genes/permcont.csv", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")

ggplot(contdata, aes(x=rows, y=data)) + geom_point()

rownames(cont)
# EdgeR plots ----
jpeg("generalfigures/BCVplot_edgeR.jpg")
plotBCV(est)
dev.off()

jpeg("generalfigures/QLDISP_edgeR.jpg")
plotQLDisp(fit)
dev.off()

jpeg("generalfigures/MDS_edgeR.jpg")
plotMDS(logcpm)
dev.off()


# jpeg("generalfigur>es/SpliceDGE_edgeR.jpg")
# plotSpliceDGE(logcpm)
# dev.off()

# end ----
library(matrixStats)
stats <- txi$abundance
tpm <- transform(stats, sd=rowSds(stats), avg=rowMeans(stats), median = rowMedians(stats))

dge <- DGEList(txi$counts)
cpm <- cpm(dge, log=TRUE)
head(cpm)
data.frame(cpm)
cpmtibble <- as_tibble(cpm, rownames="geneID")
cpmpivot <- pivot_longer(cpmtibble, cols=seq(2, ncol(cpmtibble)), names_to = "samples", values_to = "expression")



ggplot(tpm) +
  aes(x = sd, y = median) +
  geom_hex()

ggplot(tpm) +
  aes(x = sd, y = avg) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex() +
  labs(y="Median", x = "Standard deviation", title="Transcripts per million (TPM)", subtitle="unfiltered, non-normalized data") +
  theme_classic() +
  theme_dark() +
  theme_bw()

ggplot(cpmpivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",geom = "point",shape = 95,size = 10,color = "black",show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",title="Log2 Counts per Million (CPM)",subtitle="unfiltered, non-normalized",caption=paste0("produced on ", Sys.time())) +
  theme_bw()

keepers <- rowSums(cpm > 1)>=3
table(rowSums(dge$counts==0)==10)


# venndiagrams ----

library(VennDiagram)


bodygenes <- names(Deseq)
bodygenes <- bodygenes[!bodygenes %in% rownames(Leg)]
leggenes <- names(Deseq)
leggenes <- leggenes[!leggenes %in% rownames(Body)]
names(Deseq) - rownames(Body)

#
bodyup <- table[table["partleg"] > 0]


# Venn Diagram 1----
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.diagram(
  x = list(bodygenes, leggenes),
  category.names = c("Body" , "Leg"),

  filename = '14_venn_diagramm.png',
  output=FALSE,
  imagetype="png",
  compression="lzw",


  lwd = 2,
  lty = 'blank',
  fill = c("red", "blue"),
)
