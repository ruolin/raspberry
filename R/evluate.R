rproc <- function(collector, truth, bayes=T){
  gene_nominal_p = c()
  for (i in 1:length(collector)){
    if (bayes) {
      condp = collector[[i]]$pvalues[grepl(",2]", names(collector[[i]]$pvalues))]
      names(condp) = collector[[i]]$gene_ids
    } else {
      condp = collector[[i]]$pvalues[,2]
      names(condp) = collector[[i]]$gene_ids[-1]
    }
    #trans_nominal_p = c(trans_nominal_p, condp)
    uniq_ids = unique(names(condp))
    offset = 0
    for (j in 1:length(uniq_ids)) {
      gene_level = condp [ which(names(condp) == uniq_ids[j])]
      idx = offset + which.min(gene_level)
      gene_nominal_p = c(gene_nominal_p, condp[idx])
      offset = offset + length(gene_level)
    }
  }
  
  gene_pvalues = rep(1, length(unique(names(gene_nominal_p))))
  names(gene_pvalues) = unique(names(gene_nominal_p))
  for(i in 1:length(gene_nominal_p)) {
    cur_gene = names(gene_nominal_p[i])
    gene_pvalues[cur_gene] = min (gene_pvalues[cur_gene], gene_nominal_p[i])
  }
  labels = as.integer(names(gene_pvalues) %in% truth)
  list(scores = gene_pvalues, labels =labels)
}


evaluate <- function(raspberry_result, truth_file, ...) {
  truth = scan(truth_file, "")
  raspberry <- rproc(raspberry_result, truth, bayes=T)
  rroc <- pROC::plot.roc(raspberry$label, raspberry$score, main="", col= "lightslateblue", legacy.axes=TRUE,
                       xlab="False positive rate (FPR)", ylab="True positive rate (TPR)",identity=F,cex.axis=2.2, cex.lab=2.2, lwd=4)
  print(paste("raspberry AUC:", rroc$auc))
  other_results <- list(...)
  nprog <- length(other_results)
  cols <- RColorBrewer::brewer.pal(nprog, "Set1")
  
  legend_txt = c("Raspberry") 
  legend_cols = c("lightslateblue")
  for (i in 1:nprog) {
    program = names(other_results)[i]
    roc_df = prepare_roc(other_results[[i]], truth, program)
    roc = pROC::plot.roc(roc_df$label, roc_df$score, main="", col= cols[i], legacy.axes=TRUE,
                          xlab="False positive rate (FPR)", ylab="True positive rate (TPR)",identity=F,cex.axis=2.2, cex.lab=2.2, lwd=4, add=T)
    print (paste(program, "AUC:", roc$auc))
    legend_txt = c(legend_txt, program)
    legend_cols = c(legend_cols, cols[i])
  }
  
  legend("bottomright", legend = legend_txt,
         col = legend_cols, lwd = 5)
}

prepare_roc <-function(file, truth, program = "Cuffdiff2") {
  df = read.table(file, header = T)
  if(program == "Cuffdiff2") {
    df_ok = df[df$status == "OK",]
    df_ok[,"label"] = 0
    df_ok[df_ok$gene %in% truth, "label"] = 1
    list(score = df_ok$p_value, label = df_ok$label)
  }
  else if(program == "DEXSeq") {
    genes_min_p = tapply(df$pvalue, df$groupID, min)
    genes_min_p_nona = genes_min_p[ !is.na(genes_min_p) ]
    list(score = genes_min_p_nona, label = as.integer(names(genes_min_p_nona) %in% truth))
  }
  else if(program == "DSGseq") {
    df_ok = df[df$is_filtered == "FALSE",]
    df_ok[, "label"] = 0
    df_ok[df_ok$gene_name %in% truth, "label"] = 1
    list(score = df_ok$NB_stat, label = df_ok$label)
  }
  else {
    message("Unknown tool")
    stopifnot(FALSE)
  }
}


