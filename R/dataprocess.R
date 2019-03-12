getXy <- function(data, cases, include_geneid = FALSE){
  gene_idx = "gene_id"
  gene_count = "gene_frag_count"
  iso_idx = "transcripts"
  cond_prop_idx = "conditional_probabilities"
  class_prop_idx ="class_probabilities"
  path_idx = "path_symbol"
  bincount_idx = "path_count"
  fpkm_idx = "FPKMs"
  sample_idx = "sample"
  niso = 0
  if (grepl(",", data[1, class_prop_idx])) {
    niso = length(unlist(strsplit(data[1, class_prop_idx], ",")))
  } else {
    niso = 1
  }
  iso_names = NULL
  if (niso == 1) {
    iso_names = data[1, iso_idx]
  } else {
    iso_names  = unlist(strsplit(data[1, iso_idx], ","))
  }
  samples = unique(data[, sample_idx])
  fpkm_mat = matrix(nrow=length(samples), ncol = niso, NA)
  rownames(fpkm_mat) = samples
  colnames(fpkm_mat) = iso_names
  nrows = nrow(data)
  for (i in 1:nrows) {
    if ( all(is.na(fpkm_mat[data[i, sample_idx],])) ) {
      fpkm_mat[data[i, sample_idx],] = as.numeric(unlist(strsplit(data[i, fpkm_idx], ",")))
    }
  }
  cond_prop_mat = matrix(nrow=nrows, ncol=niso, 0)
  #class_prop_mat = matrix(nrow=nrows, ncol=niso, 0)
  for (i in 1:nrows) {
    #if ( length(unlist(strsplit(data[i,cond_prop_idx], ","))) != niso) { 
      # Right now cannot handle isoforms that do not overlap
      # For example gene ARHGEF7
      #return (list(data=data, niso = 1, nbins = 1))
    #}
    if (niso == 1) {
      cond_prop_mat[i,] = as.numeric(data[i, cond_prop_idx])
    } else {
      cond_prop_mat[i,] = as.numeric(unlist(strsplit(data[i, cond_prop_idx], ",")))
    }
  }
  #Filter row contains all 0 thetas
  row_to_remove = c()
  for (i in 1:nrows) {
    if (sum(cond_prop_mat[i,]) < 1e-10) {
      row_to_remove = c(row_to_remove, i)
    }
  }
  
  if (length(row_to_remove) != 0) {
    cond_prop_mat = cond_prop_mat[-row_to_remove]
    #class_prop_mat = class_prop_mat[-row_to_remove]
    data = data[-row_to_remove,]
    nrows = nrow(data)
  }
  
  nbins = length(unique(data[,path_idx]))
  #if (niso == 1  || nbins == 1) {
  if (nbins == 1) {
    return (list(data=data, niso = niso, nbins = nbins))
  }
  
  #print(head(data))
  # sample matrix
  sm = model.matrix( ~sample-1, data=data) 
  nsample = ncol(sm)
  sm = sm[,-ncol(sm)]
  X = sm
  
  # iso matrix
  iso_matrix = matrix(0, nrow = nrows, ncol=niso) 
  colnames(iso_matrix) = iso_names
  X = cbind(X, iso_matrix) 
  X = cbind(X, model.matrix(~data[, path_idx] -1)[,-1]) 
  
  # update for condition
  for (k in 1:nrows) {
    if (data[k,1] %in% cases) {
      sm[k,1] = 1
    } else {
      sm[k,1] = 0
    }
  }
  colnames(sm)[1] = "condition"
  ###
  ### expand the model matrix for poisson regression
  ### contain indicator variable for each isoform
  ###
  iso_start_idx = which(iso_idx == colnames(data))
  large_X = X
  large_sm = sm
  #stopifnot(niso > 1)
  
  for (i in 2:niso) {
    large_X = rbind(large_X,X)
    large_sm = rbind(large_sm, sm)
  }
  
  
  for (i in 1:nrows) {
    for (j in 1:niso) {
      idx = i + (j -1)*nrow(data)
      large_X[idx, iso_start_idx + j] = 1
    }
  }
  large_X = large_X[,-(iso_start_idx+1)] # identifiability constraints
  
  ### the interaction effects
  for (i in 1:(niso-1)) {
    #colnames(large_X)[iso_start_idx+i]
    tmp_X=large_X[,iso_start_idx+i] * large_sm
    colnames(tmp_X)[1] = "condition"
    colnames(tmp_X) = paste0(colnames(large_X)[iso_start_idx+i], "*", colnames(tmp_X))
    large_X = cbind(tmp_X, large_X)    
    iso_start_idx = iso_start_idx + ncol(tmp_X)
  }
  
  ## SX is for multinomial regression
  row.names(large_X) = rep(1:nrow(large_X))
  list(X=large_X, data=data, SX=sm, niso = niso, nbins = nbins, nsample=nsample, fpkm_mat = fpkm_mat)
}
#############
#############
#############

GetPiMatrix <-function(df_, RANDOM_PI_START) {
  pi = NULL
  g2t_map = list()
  t2g_map = list()
  class_prop_idx = "class_probabilities"
  for (i in 1:length(df_[[1]])) {
    gene_pi = NULL
    gene_name = names(df_[[1]])[i]
    dat = df_[[1]][[i]]
    niso = 0
    if (grepl(",", dat[1, class_prop_idx])) {
      niso = length(unlist(strsplit(dat[1, class_prop_idx], ",")))
    } else {
      niso = 1
    }
    
    if (RANDOM_PI_START) {
      gene_pi = array(1/niso, dim = c(1, niso))
    } else {
      gene_pi = array(as.numeric(unlist(strsplit(dat[1, class_prop_idx], ","))), dim = c(1, niso))
    }
    rownames(gene_pi) = dat[1,"sample"]
    if (niso == 1) {
      colnames(gene_pi) = dat[1,"transcripts"]
      g2t_map[[length(g2t_map) + 1]] = dat[1,"transcripts"]
      names(g2t_map)[length(g2t_map)] = gene_name
      t2g_map[[length(t2g_map) + 1]] = gene_name
      names(t2g_map)[length(t2g_map)] = dat[1,"transcripts"]
    } else {
      colnames(gene_pi) = unlist(strsplit(dat[1,"transcripts"], ","))
      txs = unlist(strsplit(dat[1,"transcripts"], ","))
      g2t_map[[length(g2t_map) + 1]]  =  txs
      names(g2t_map)[length(g2t_map)] = gene_name
      for (tx in txs) {
        t2g_map[[length(t2g_map) + 1]] = gene_name
        names(t2g_map)[length(t2g_map)] = tx
      }
    }
    skip_gene = FALSE
    for (k in 2:length(df_)) {
      if (gene_name %in% names(df_[[k]])) {
        other_dat = df_[[k]][[gene_name]]
        if (RANDOM_PI_START) {
          gene_pi = rbind(gene_pi, matrix(1/niso, ncol = niso, nrow = 1))
        } else {
          if (niso == 1) {
            gene_pi = rbind(gene_pi, as.numeric(other_dat[1, class_prop_idx]))
          } else {
            gene_pi = rbind(gene_pi, as.numeric(unlist(strsplit(other_dat[1, class_prop_idx], ","))))
          }
        }
        rownames(gene_pi)[k] = other_dat[1,"sample"]
      }
      else {
        skip_gene = TRUE
        break
      }
    }
    if (skip_gene != TRUE) {
      if (is.null(pi)) {
        pi = gene_pi
      } else {
        pi = cbind(pi, gene_pi)
      }
    }
  }
  list(pi=pi, g2t_map = g2t_map, t2g_map = t2g_map)
}


###########
###########

SubsetGeneIds <- function(genelist, e, s=1){
  if (s > length(genelist)) return (NULL)
  if (e > length(genelist)) e = length(genelist)
  gids = genelist[s:e]
  gids
}

GetAllGeneIds <- function(file) {
  df = read.table(file, header=T, stringsAsFactors = F)
  ids = unique(df$gene_id)
}


FilterBySingleSampleEst <- function(pi, cases, CUTOFF = 0.3){
  result = c()
  colidx = rownames(pi$pi) %in% cases
  colnames(pi$pi)[i]
  for( i in 1:dim(pi$pi)[2]) {
    group1 = pi$pi[colidx,i]
    group2 = pi$pi[!colidx,i]
    if (mean(group1) == mean(group2)) {
      next
    }
    ttest = t.test(group1, group2)
    if ( ttest$p.value < CUTOFF) {
      result = c(result, pi$t2g_map[colnames(pi$pi)[i]])
    }
  }
  unique(unlist(result))
}
