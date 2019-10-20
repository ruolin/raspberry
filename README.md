# Raspberry
RNA-Seq differential spliced transcripts toolkit

Raspberry utilizes a multinomial logistic regression to regress transcript relative abundances on a set of covariates to explain replicate and condition effects.
To better overcome the coverage bias, Raspberry uses a dual-phase algorithm.
During the bias correction phase, another multinomial logistic regression model is trained on single-isoform loci across all samples to discover relationships between the subexon path probability and the underlying local sequence information such as GC-content, existence of high GC-stretches and hexamer context. 
Then during the second phase, the differential splicing analysis algorithm uses the fitted path probabilities instead of the observed path probabilities. To account for the count overdispersion, Raspberry employs an empirical Bayesian model which places shrinkage priors on the multinomial logistic regression coefficients, which represent isoform-level dispersion and etc..
These priors are estimated by borrowing information across all loci and all samples.

To our knowledge, Raspberry is the first to simultaneously estimate transcript abundances and identify differential alternative splicing at the transcript level.
Existing methods either measure at the level of splicing events (for example inclusion or exclusion of a particular cassette exon), e.g., DEXSeq,  or detect at the level of genes, e.g., Cuffdiff 2.
In addition, Raspberry combines a bias correction step into the detection of differential alternative splicing.
Although many transcript quantification methods employ bias correction steps, the effect of coverage bias correction has not been observed and fully studied in the context of differential alternative splicing detection until Raspberry.

### Detection of differential alternative splicing from two groups of RNA-SEQ samples.

prerequisites
================================

1. Strawberry. This can be found at https://github.com/ruolin/strawberry


2. devtools.
```
install.packages("devtools")
```

3. Two groups of BAM files.

4. A genome file.

Preparing the input for rasphberry
======================================
Strawberry is used to preprocess the bam files and calculate the exon_bin counts. 
run strawberry on each individual bams with a given gff file and genome file and `-f` option. 
  
    strawberry ERR188021.sorted.bam -g /path/to/your/annotation.gff -b /path/to/your/genome.fa -o ERR188021.gtf -f frag_context_ERR188021.txt
  
For non-model organism, can run Strawberry to reconstruct gene models and then merged them by Cuffmerge. 
  
    1. strawberry ERR188021.sorted.bam -o ERR188021.gtf
  
    2. strawberry ERR188088.sorted.bam -o ERR188088.gtf
  
    3. echo ERR188021.gtf > samples.lst
  
    4. echo ERR188088.gtf >> samples.lst
  
    5. Cuffmerge samples.lst

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
