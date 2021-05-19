txi <- function(directory, quantdatafoldername, metadata, fishtype){
  library(tximport)
  setwd(directory)
  metadata <- read.csv(metadata, header=TRUE)
  files <- file.path(directory, quantdatafoldername, metadata$name)
  names(files) <- metadata$name
  genes <- read.table("genes2.tabular")
  txi <- tximport(files, type=fishtype, tx2gene=genes)
  return(list(txi, metadata))
}

initialDeseq <- function(txi, metadata, design=~part + pesticides){
  library(DESeq2)
  # Deseq2 dataset from tximport
  data <- DESeqDataSetFromTximport(txi, colData=metadata, design = design)
  Deseq <- DESeq(data)
  results <- results(Deseq)
  print(summary(results))
  return(list(Deseq, results))
}

reDeseq <- function(dds){
  ddsMF <- dds
  design(ddsMF) <- formula(~pesticides)
  ddsMF <- DESeq(ddsMF)
  Mres <- results(ddsMF)
  return(list(ddsMf, Mres))
  }

initialEdgeR <- function(txi, meta){
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
  
  pest <- factor(meta$pesticides, levels=c("perm", "deet", "control"))
  part <- factor(meta$part, levels=c("leg", "body"))
  design <- model.matrix(~0 + pest + part)
  keep <- filterByExpr(y, design)
  EdgeR <- y[keep, ]
  return(list(EdgeR, design))
}

estEdgeR <- function(EdgeR, design){
  library(edgeR)
  y <- estimateGLMCommonDisp(EdgeR, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmQLFit(y, design) 
  
  return(fit)
}
