getHiddenTable <- function(obs_table, pi, theta) {
  niso = 0
  nsample = 0
  if (is.vector(pi)) {
    niso = 1
    nsample = length(pi)
  } else {
    niso = ncol(pi)
    nsample = nrow(pi)
  }
  npath = ncol(obs_table)
  hidden_table = array(0, dim=c(nsample, npath, niso))
  for (j in 1:nsample) {
    for (l in 1:npath) {
      if (niso == 1) {
        hidden_table[j, l, 1] = obs_table[j, l]
      } else {
        den = theta[l, ] %*% pi[j,]
        for (k in 1:niso) {
          if (den == 0) {
            hidden_table[j, l, k] = 0
          } else {
            hidden_table[j, l, k] = obs_table[j, l]  * pi[j, k] * theta[l, k] / den
          }
        }
      }
    }
  }
  hidden_table
}


getTheta <- function(hidden_table) {
  NUMERIC_LOWEST = 1e-14
  npath = dim(hidden_table)[2]
  niso = dim(hidden_table)[3]
  #nsample = dim(hidden_table)[1]
  theta = array(0, dim=c(npath, niso))
  for (k in 1:niso) {
    denom = sum(hidden_table[,,k])
    if (denom < NUMERIC_LOWEST) {
      theta[,k] = 0
    }else {
      for (l in 1:npath) {
        theta[l,k] = sum(hidden_table[,l,k])/denom
      }
    }
  }
  theta
}

getY <- function(hidden_table) {
  #npath = dim(hidden_table)[2]
  niso = dim(hidden_table)[3]
  nsample = dim(hidden_table)[1]
  Y = matrix(0, nrow = nsample, ncol = niso)
  for (j in 1:nsample) {
    for (k in 1:niso) {
      Y[j, k] = sum(hidden_table[j,,k])
    }
  }
  Y
}

##########

GenewiseHiddenTable <- function (mat, biasfit, pi, normalize) {
  ### initialize
  path_idx = "path_symbol"
  sample_idx = "sample"
  count_idx = "path_count"
  gc_idx = "path_gc_content"
  entropy_idx = "path_hexmer_entropy"
  iso_idx = "transcripts"
  theta_idx = "conditional_probabilities"

  paths = unique(mat$data[,path_idx])
  samples = unique(mat$data[,sample_idx])
  npaths = length(paths)
  nsample = mat$nsample
  niso = mat$niso

  ## Observed table
  obs_table = array(0, dim=c(nsample, npaths))
  colnames(obs_table) = paths
  rownames(obs_table) = samples
  for ( i in 1:nrow(mat$data)) {
    obs_table[mat$data[i,sample_idx], mat$data[i,path_idx]] = as.numeric(mat$data[i,count_idx])
  }
  # Calculate pi
  stopifnot(!is.null(pi))

  #calculate theta
  #everage theta across all samples
  theta = array(0, dim=c(npaths, niso))
  rownames(theta) = paths
  path_counts = table(mat$data[path_idx])
  for (i in 1:nrow(mat$data)) {
    if (niso == 1) {
      theta[mat$data[i, path_idx], ] =  theta[mat$data[i, path_idx], ] + as.numeric(mat$data[i, theta_idx])
    } else {
      theta[mat$data[i, path_idx], ] =  theta[mat$data[i, path_idx], ] + as.numeric(unlist(strsplit(mat$data[i, theta_idx], ",")))
    }
  }

  for (i in 1:nrow(theta)) {
    theta[i,] = theta[i, ] / path_counts[rownames(theta)[i]]
  }
  if (!is.null(biasfit)) {
    fittedbias = getFittedBias(mat$data, biasfit)
    if (!is.null(fittedbias)) {
      new_theta = theta
      for(j in 1:ncol(theta)) {
        nonzeros = which(theta[,j]>0)
        if (normalize) {
          denom = sum(fittedbias$fitted[nonzeros])
        } else {
          denom = sum(fittedbias$fitted)
        }
        for (i in 1:length(fittedbias$fitted)) {
          if (i %in% nonzeros) {
            new_theta[i,j] = fittedbias$fitted[i]/denom
          } else {
            new_theta[i,j] = 0
          }
        }
      }
      theta = new_theta
    }
  }

  ##
  hidden_table = getHiddenTable(obs_table, pi, theta)
  Y = getY(hidden_table)
  X_mul = unique(mat$SX)
  SX = mat$SX
  if (niso == 1) {
    colnames(Y) = mat$data[1,iso_idx]
  } else {
    colnames(Y) = unlist(strsplit(mat$data[1,iso_idx], ","))
  }
  rownames(Y) = rownames(obs_table)
  list(y_mul=Y, X_mul = X_mul, theta = theta, obs_table = obs_table)
}


##########
##########
