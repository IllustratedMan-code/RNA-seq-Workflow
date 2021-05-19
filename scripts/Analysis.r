

setwd("/home/david/Documents/BenoitLab/RNA-seq/")
source("scripts/functions.r")
directory <- "/home/david/Documents/BenoitLab/R/largeRNASEQ/"
quantdatafoldername <- "CorrectedQuantData"
fishtype <- "sailfish"

a <- txi(directory, quantdatafoldername, "metadata.csv", fishtype)
#De <- initialDeseq(a[[1]], a[[2]])
Ed <- initialEdgeR(a[[1]], a[[2]])
Eest <- estEdgeR(Ed[[1]], Ed[[2]])
lrt <- glmLRT(Eest, coef=1:4)
test <- glmQLFTest(Eest, coef=1:4)

table <- test$table
body <- table[table$logFC.partbody > 2, ]
leg <- table[table$logFC.partbody < -2, ]

setwd("/home/david/Documents/BenoitLab/RNA-seq/")
write.table(table, file="edgeR-Genes/TotalExpr.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(leg, file="edgeR-Genes/legExpr.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")
write.table(body, file="edgeR-Genes/bodyExpr.csv", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")

deetvcontrol <- table$logFC.pestdeet - table$logFC.pestcontrol
permvcontrol <- table$logFC.pestperm - table$logFC.pestcontrol
deetvcontrolvbody <- body$logFC.pestdeet - body$logFC.pestcontrol
permvcontrolvbody <- body$logFC.pestperm - body$logFC.pestcontrol
deetvcontrolvleg <- leg$logFC.pestdeet - leg$logFC.pestcontrol
permvcontrolvleg <- leg$logFC.pestperm - leg$logFC.pestcontrol

table$deetvcontrol <- deetvcontrol
table$permvcontrol <- permvcontrol
body$deetvcontrol <- deetvcontrolvbody
body$permvcontrol <- permvcontrolvbody
leg$deetvcontrol <- deetvcontrolvleg
leg$permvcontrol <- permvcontrolvleg

updeet <- table[table$deetvcontrol > 2,]
downdeet <- table[table$deetvcontrol < -2,]



updeetbody <- body[body$deetvcontrol > 2,]
downdeetbody <- body[body$deetvcontrol < -2,]
updeetleg <- leg[leg$deetvcontrol > 2,]
downdeetleg <- leg[leg$deetvcontrol < -2,]


upperm <- table[table$permvcontrol > 2,]
downperm <- table[table$permvcontrol < -2,]
uppermbody <- body[body$permvcontrol > 2,]
downpermbody <- body[body$permvcontrol < -2,]
uppermleg <- leg[leg$permvcontrol > 2,]
downpermleg <- leg[leg$permvcontrol < -2,]

upperm <- table[table$permvcontrol > 2,]
downperm <- table[table$permvcontrol < -2,]
uppermbody <- body[body$permvcontrol > 2,]
downpermbody <- body[body$permvcontrol < -2,]
uppermleg <- leg[leg$permvcontrol > 2,]
downpermleg <- leg[leg$permvcontrol < -2,]

uppermgenes <- rownames(upperm)
downpermgenes <- rownames(downperm)
uppermleggenes <- rownames(uppermleg)
downpermleggenes <- rownames(downpermleg)
uppermbodygenes <- rownames(uppermbody)
downpermbodygenes <- rownames(downpermbody)

updeetgenes <- rownames(updeet)
downdeetgenes <- rownames(downdeet)
updeetleggenes <- rownames(updeetleg)
downdeetleggenes <- rownames(downdeetleg)
updeetbodygenes <- rownames(updeetbody)
downdeetbodygenes <- rownames(downdeetbody)

genelist <- list(upPerm = uppermgenes, downPerm = downpermgenes,
                upPermBody=uppermbodygenes, downPermBody=downpermbodygenes,
                upPermLeg = uppermleggenes, downPermLeg = downpermleggenes,
                upDeet = updeetgenes, downDeet = downdeetgenes,
                upDeetLeg = updeetleggenes, downDeetLeg = downdeetleggenes,
                upDeetBody = updeetbodygenes, downDeetBody = downdeetbodygenes, table=rownames(table), leg=rownames(leg), body=rownames(body))

for (i in seq(length(genelist))){
    namelist <- names(genelist)
    write.table(genelist[i], file=paste("Genes/Genelists","/", namelist[i], ".csv", sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
}
