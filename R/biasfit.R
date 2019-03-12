getBiasMat <- function(gene_df){
  gene_idx = "gene_id"
  gene_count_idx = "gene_frag_count"
  iso_idx = "transcripts"
  exonbin_idx = "path_symbol"
  bincount_idx = which(colnames(gene_df) == "path_count")
  
  iso_names  = unlist(strsplit(gene_df[,iso_idx], ","))
  niso = length(iso_names)
  nrows = nrow(gene_df)
  nbins = length(unique(gene_df[,exonbin_idx]))
  bias_mat = gene_df[,bincount_idx:ncol(gene_df)]
  bias_mat = cbind(gene_df[, c(gene_idx, gene_count_idx)], bias_mat)
  return (list(bias_mat = bias_mat, niso=niso, nbins = nbins))
}

getBiasCoef<- function(total_bias_mat_){
    head(total_bias_mat_)
   #total_bias_mat_ = head(total_bias_mat, n=200)
   require(splines)
   gc_knots = c(0.4,0.5,0.6)
   gc_boundary_knots = c(0.3, 0.7)
   entropy_knots = c(4,5,6)
   entropy_boundary_knots = c(3,7)

   gc_stretch_idx = 6:9
   gene_mat = model.matrix(~total_bias_mat_[,1]-1)
   #colnames(total_bias_mat_)[1:9] = c("gene","gene_count","bin_count", "gc_content", "entropy", "high_gc_20_0.8","high_gc_20_0.9","high_gc_40_0.8","high_gc_40_0.9") 
   gc_cubic = ns (total_bias_mat_[, "path_gc_content"], knots = gc_knots, Boundary.knots = gc_boundary_knots)
   entropy_cubic = ns(total_bias_mat_[, "path_hexmer_entropy"], knots = entropy_knots, Boundary.knots = entropy_boundary_knots)
   colnames(gc_cubic) = paste0("gc_cubic_basis", 1:ncol(gc_cubic))
   colnames(entropy_cubic) = paste0("entropy_cubic_basis", 1:ncol(entropy_cubic))
   mat = cbind(total_bias_mat_[,"path_count"], gc_cubic, entropy_cubic, total_bias_mat_[, gc_stretch_idx], gene_mat)
   colnames(mat)[1] = "bin_count"
   pr.fit = glm(bin_count ~ . - bin_count - 1, data = mat, family = "poisson")
   num_coef = ncol(gc_cubic) + ncol(entropy_cubic) + length(gc_stretch_idx)
   list(coef = coef(pr.fit)[1:num_coef], num_gc_basis = ncol(gc_cubic), num_entropy = ncol(entropy_cubic), fit = pr.fit, covar= mat[,2:9], bin_count = total_bias_mat_[,3])
   #list(coef = coef(pr.fit)[1:num_coef], num_gc_basis = ncol(gc_cubic), num_entropy = ncol(entropy_cubic))
}

getFittedBias <- function(data, biasfit){
  #data=res$data
  logy = NULL
  path_idx = "path_symbol"
  gc_idx = "path_gc_content"
  entropy_idx = "path_hexmer_entropy"
  gc_stretch_idx = c("gc_strecth_0.8_20", "gc_strecth_0.9_20", "gc_strecth_0.8_40", "gc_strecth_0.9_40")
  gc_knots = c(0.4,0.5,0.6)
  gc_boundary_knots = c(0.3, 0.7)
  entropy_knots = c(4,5,6)
  entropy_boundary_knots = c(3,7)
  gc_cubic = ns (data[,gc_idx], knots = gc_knots, Boundary.knots = gc_boundary_knots)
  entropy_cubic = ns (data[,entropy_idx], knots = entropy_knots, Boundary.knots = entropy_boundary_knots)
  colnames(gc_cubic) = paste0("gc_cubic_basis", 1:ncol(gc_cubic))
  colnames(entropy_cubic) = paste0("entropy_cubic_basis", 1:ncol(entropy_cubic))
  bias_covar_dup = cbind(gc_cubic, entropy_cubic, data[gc_stretch_idx])
  bias_covar_dup[,ncol(bias_covar_dup)+1] = data[,path_idx]
  bias_covar = bias_covar_dup[!duplicated(bias_covar_dup[,ncol(bias_covar_dup)]),]
  rownames(bias_covar) = bias_covar[,ncol(bias_covar)]
  bias_covar = bias_covar[,-ncol(bias_covar)]
  na_idx = which(is.na(biasfit$coef))
  if (length(na_idx) == 0) {
    logy = as.matrix(bias_covar)%*% biasfit$coef  
  } else {
    logy = as.matrix(bias_covar[,-na_idx])%*% biasfit$coef[-na_idx]
  }
  list(fitted = exp(logy))
}


