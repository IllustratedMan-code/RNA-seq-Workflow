library("tidyverse")
library("ggpubr")

EdgeR <- read.csv("/home/david/Documents/BenoitLab/RNA-seq/edgeR-Genes/bodyExpr.csv")
Deseq <- read.csv("/home/david/Documents/BenoitLab/RNA-seq/Deseq-Genes/total.csv")
common <- intersect(row.names(Deseq), row.names(EdgeR))
EdgeRCommon <- EdgeR[common, ]
DeseqCommon <- Deseq[common, ]

ECol <- EdgeRCommon["logFC.partbody"]
DCol <- -1*DeseqCommon["log2FoldChange"]
newD <- data.frame(DCol, ECol)

library("ggpubr")
pearson <- ggscatter(newD, x = "log2FoldChange", y = "logFC.partbody",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Deseq", ylab = "EdgeR")
pearson
ggsave("/home/david/Documents/BenoitLab/RNA-seq/figures/figure1/pearson.png", pearson, png())
ggsave("/home/david/Documents/BenoitLab/RNA-seq/figures/figure1/pearson.pdf", pearson, pdf())

library(VennDiagram)
ECol <- row.names(EdgeR["logFC.partbody"])
DCol <- row.names(-1*Deseq["log2FoldChange"])
venn.diagram(
  x = list(na.omit(ECol), na.omit(DCol)),
  category.names = c("EdgeR", "Deseq"),
  filename = "/home/david/Documents/BenoitLab/RNA-seq/figures/figure1/vennBody.png",
  output=TRUE
)
