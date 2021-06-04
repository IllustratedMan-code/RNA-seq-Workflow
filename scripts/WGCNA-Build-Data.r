setwd("/home/david/Documents/BenoitLab/R/largeRNASEQ/CorrectedQuantData/")
quantfiles <- list.files(".", pattern = ".tabular")
df <- data.frame(row.names = read.delim(quantfiles[1])$Name)
for (fileName in quantfiles) {
  file <- read.delim(fileName)
  df[tools::file_path_sans_ext(fileName)] <- file$TPM
}
write.csv(df, "../WGCNACorrectedQuant.csv", quote = FALSE)
head(df)

setwd("/home/david/Documents/BenoitLab/R/largeRNASEQ")
metadata <- read.csv("metadata.csv", row.names = 1)
truthTable <- data.frame(row.names = row.names(metadata))
for (column in colnames(metadata)) {
  for (row in unique(metadata[column])) {
    for (item in row) {
      truthTable[, paste(column, item, sep = "_")] <- as.integer(item == metadata[[column]])
    }
  }
}
truthTable <- data.frame(sample = row.names(truthTable), truthTable)
write.csv(truthTable, "WGCNATruthTable.csv", row.names=FALSE)
truthTable

setwd("/home/david/Documents/BenoitLab/R/largeRNASEQ")
truthTable <- data.frame(sample = colnames(df))
truthTable$Control <- c(
  1, 1, 1, 1, 1, 1,
  1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0
)
truthTable$DEET <- c(
  0, 0, 0, 0, 0, 0,
  0, 1, 1, 1, 1, 1,
  1, 0, 0, 0, 0, 0
)
truthTable$Permethrin <- c(
  0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0,
  0, 1, 1, 1, 1, 1
)
truthTable$Body <- c(
  1, 1, 1, 0, 0, 0,
  0, 1, 1, 1, 0, 0,
  0, 1, 1, 0, 0, 0
)
truthTable$Leg <- c(
  0, 0, 0, 1, 1, 1,
  1, 0, 0, 0, 1, 1,
  1, 0, 0, 1, 1, 1

write.csv(truthTable, "./WGCNATruthTable.csv", row.names = FALSE)

A data.frame: 6 × 18

| <!--/--> | Cont_Bod_1 &lt;dbl&gt; | Cont_Bod_2 &lt;dbl&gt; | Cont_Bod_3 &lt;dbl&gt; | Cont_Leg_1 &lt;dbl&gt; | Cont_Leg_2 &lt;dbl&gt; | Cont_Leg_2a &lt;dbl&gt; | Cont_Leg_3 &lt;dbl&gt; | Deet_Bod_1 &lt;dbl&gt; | Deet_Bod_2 &lt;dbl&gt; | Deet_Bod_3 &lt;dbl&gt; | Deet_Leg_1 &lt;dbl&gt; | Deet_Leg_2 &lt;dbl&gt; | Deet_Leg_3 &lt;dbl&gt; | Per_Bod_1 &lt;dbl&gt; | Per_Bod_3 &lt;dbl&gt; | Per_Leg_1 &lt;dbl&gt; | Per_Leg_2 &lt;dbl&gt; | Per_Leg_3 &lt;dbl&gt; |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dv.1 | 22.2763 | 37.8193 | 40.3742 | 61.417900 | 0.0681662 | 83.278200 | 69.070000 | 48.1591 |  43.6397 | 50.0788 | 86.710000 | 94.912900 | 85.979000 | 52.7464 | 56.6623 | 9.02202e+01 | 96.067400 | 76.999000 |
| dv.2 | 15.9128 | 44.7599 | 54.8798 |  0.786553 | 0.0756211 |  0.722977 |  0.493229 | 93.6505 | 105.8530 | 72.7003 |  2.753670 |  0.870014 |  3.473720 | 73.2376 | 80.6752 | 1.26751e+00 |  0.855653 |  0.255794 |
| dv.3 | 10.6041 | 11.7978 | 14.3419 |  9.187200 | 0.0756683 | 14.169700 | 12.297800 | 16.7729 |  14.5622 | 18.3525 | 14.070000 | 17.049100 | 15.897900 | 16.5017 | 18.0218 | 1.62399e+01 | 18.330600 | 14.676500 |
| dv.4 | 16.5992 | 19.8922 | 23.9323 | 24.027200 | 0.0794267 | 30.721800 | 27.661300 | 24.9121 |  21.8748 | 27.2285 | 29.715900 | 35.491900 | 34.275800 | 23.6021 | 25.1121 | 2.79280e+01 | 31.287200 | 25.525300 |
| dv.5 |  0.0000 |  0.0000 |  0.0000 |  0.294700 | 0.0834384 |  0.875706 |  0.715919 |  0.0000 |   0.0000 |  0.0000 |  0.320064 |  0.391045 |  0.193081 |  0.0000 |  0.0000 | 2.79643e-05 |  0.671768 |  0.217907 |
| dv.6 |  0.0000 |  0.0000 |  0.0000 |  0.484270 | 0.0837346 |  0.431047 |  0.759755 |  0.0000 |   0.0000 |  0.0000 |  0.329439 |  0.940694 |  1.231160 |  0.0000 |  0.0000 | 1.37827e+00 |  0.962562 |  0.141731 |
