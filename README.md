# raspberry
RNA-Seq differential spliced transcripts toolkit

#prerequisites
devtools.
```
install.packages("devtools")
```

#Installation
```
library(devtools)
install_github("ruolin/raspberry")
```

#Examples
```
library(raspberry)
?ad_rd100
cases = c("/home/ruolin/Research/strawberry_comp/RD100.high.dm/hisat2/hisat_1.sorted.bam", "/home/ruolin/Research/strawberry_comp/RD100.high.dm/hisat2/hisat_2.sorted.bam", "/home/ruolin/Research/strawberry_comp/RD100.high.dm/hisat2/hisat_3.sorted.bam")
diffiso(ad_rd100, cases, FPKM_filter=1)
```
