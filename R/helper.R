
print_pred_cov<- function(gene_rf, bias_coef, gname = "NUP107") {
  target = getFittedBias(gene_rf[[1]]$gname, bias_coef)
  biasResult = cbind(rownames(target$fitted), target$fitted/(sum(target$fitted)), rep(gname, dim(target$fitted)[1]))
  colnames(biasResult) = c("path_symbol", "path_count", "gene_id")
  write.table(file =gname, biasResult, quote = FALSE, sep="\t")
}

getPvalues <-function(collector, verbose = F, bayes = T) {
  pvalues = c()
  for (i in 1:length(collector)) {
    if (bayes) {
      condp = collector[[i]]$pvalues[grepl(",2]", names(collector[[i]]$pvalues))]
      names(condp) = collector[[i]]$gene_ids
    } else {
      condp = collector[[i]]$pvalues[,"Xcondition"]
      names(condp) = rep(collector[[i]]$gene_ids[[1]], length(condp))
    }
    pvalues = c(pvalues, condp)
    if (verbose) {
      print("called")
      print(unique(dsgs)[!is.na(unique(dsgs))])
    }
  }
  list(pvalues=pvalues)
}

F1summary<-function(result, truth, critical = 0.01, verbose=F, bayes = T) {
  called = 0
  FPs = c()
  TPs = c()
  Positives = c()
  for (i in 1:length(result)) {
    dsgs = NULL
    if (bayes) {
      condp = result[[i]]$pvalues[grepl(",2]", names(result[[i]]$pvalues))]
      names(condp) = result[[i]]$gene_ids
      if (verbose) {
        print(condp)
      }
      dsgs = result[[i]]$gene_ids[condp < critical]
    } else {
      dsgs = result[[i]]$gene_ids[result[[i]]$pvalues[,"Xcondition"] < critical]
    }
    if (verbose) {
      print("called")
      print(unique(dsgs)[!is.na(unique(dsgs))])
    }

    called = called + length(unique(dsgs)[!is.na(unique(dsgs))])
    tp = intersect(truth, dsgs)

    TPs = c(TPs, tp)
    FPs = c(FPs, setdiff(dsgs, truth))
    if (verbose) {
      print("truth")
      print(intersect(truth, result[[i]]$gene_ids))
    }
    pos = intersect(truth, result[[i]]$gene_ids)
    Positives = c(Positives, pos)
  }

  a = length(TPs)/length(Positives)
  b = length(TPs)/called
  F1 = 2*a*b/(a+b)
  summary = list(TP = length(TPs), P = length(Positives), Called = called, F1 = F1)
  list(summary = summary, Recall = a, Precision = b, FPs=FPs, TPs = TPs, FNs = setdiff(Positives, TPs))
}

#########
#########

GetInitTheta<-function(df, bias_coef) {
  #df = gene_rf[[1]]$USF2
  fittedbias = getFittedBias(df, bias_coef)
  prob_str = df$conditonal_probabilites
  niso = length(unlist(strsplit(df$transcripts[1],",")))
  theta = matrix(0, ncol=niso, nrow = length(prob_str))
  for (i in 1:length(prob_str)) {
    theta[i,] = as.numeric(unlist(strsplit(prob_str[i], ",")))
  }
  new_theta = theta
  for(j in 1:ncol(theta)) {
    nonzeros = which(theta[,j]>0)
    denom = sum(fittedbias$fitted[nonzeros])
    for (i in 1:length(fittedbias$fitted)) {
      if (i %in% nonzeros) {
        new_theta[i,j] = fittedbias$fitted[i]/denom
      } else {
        new_theta[i,j] = 0
      }
    }
  }
  list(old_theta = theta, new_theta = new_theta)
}

getAllBetas <- function(coll) {
  betas = c()
  for(i in 1:length(coll)) {
    betas = c(betas, as.vector(coll[[i]]$betas))
  }
  betas
}

