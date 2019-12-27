setClass("RGene", representation(gene_id = "character", pi="matrix", 
                                 pvalue="numeric", 
                                 adj_pvalue="numeric"))
setClass("RTable", representation(significance = "data.frame", 
                                  abundance = "data.frame"))