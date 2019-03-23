# raspberry
RNA-Seq differential spliced transcripts toolkit

## prerequisites
devtools.
```
install.packages("devtools")
```

## Installation
```
library(devtools)
install_github("ruolin/raspberry")
```

## Examples
```
library(raspberry)
data("ad_rd100")
cases = c("/home/ruolin/Research/strawberry_comp/RD100.high.dm/hisat2/hisat_1.sorted.bam", "/home/ruolin/Research/strawberry_comp/RD100.high.dm/hisat2/hisat_2.sorted.bam", "/home/ruolin/Research/strawberry_comp/RD100.high.dm/hisat2/hisat_3.sorted.bam")
result = diffiso(ad_rd100, cases, FPKM_filter = 1)

true_as_gene <- system.file("extdata", "trueASgenes.txt", package="raspberry")
dsg <- system.file("extdata", "DSGresult.txt", package="raspberry")
dexseq <- system.file("extdata", "dexseq.txt", package="raspberry")
cuffdiff <- system.file("extdata", "splicing.diff", package="raspberry")
evaluate_gene_levels(result, truth_file=true_as_gene, DSGseq = dsg, DEXSeq = dexseq, Cuffdiff2 = cuffdiff)
```
Check the PDF file return by the above function.
