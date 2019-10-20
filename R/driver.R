preprocess <-function(samples_, MAX_ROW = 5000) {
  total_bias_mat = NULL
  for (i in 1:length(samples_)) {
    for (j in 1:length(samples_[[i]])) {
      mat = getBiasMat(samples_[[i]][[j]])
      if (mat$niso == 1 || mat$nbins == 1) {
        gene_name = names(samples_[[i]])[j]
        total_bias_mat = rbind(total_bias_mat, mat$bias_mat)
        if (nrow(total_bias_mat) > MAX_ROW) {
          break
        }
      }
    }
  }
  total_bias_mat
}

prepareData <- function(samples_, cases, initpi, BAYES, bias_correcter_ = NULL, target_genes = NULL, niso_min = 2, FPKM_filter = 1.0, NORMALIZE = F, verbose = F) {
  # if num_genes = k, break for only run analyses on first k genes
  # k = 0 means run for all genes
  Y = NULL
  X = NULL
  obs_table = NULL
  pi = NULL
  theta = NULL
  gene_ids = NULL
  tx_idx = "transcripts"
  for (gid in target_genes) {
    dat = NULL
    skip_gene = FALSE
    for (k in 1:length(samples_)) {
      if (is.null(samples_[[k]][[gid]])) {
        skip_gene = TRUE
        if (verbose) {
          message(paste("gene", gid, "does not express in all samples."))
        }
        break
      }
      fpkms = as.numeric(unlist(strsplit(samples_[[1]][[gid]][1, "FPKMs"], ",")))
      frac = fpkms / sum(fpkms)
      if (is.null(dat)){
        dat = samples_[[k]][[gid]]

      } else {
        dat = rbind(dat, samples_[[k]][[gid]])
      }
    }
    if (skip_gene != TRUE) {
      mat = getXy(dat, cases)
      if (mat$nbins == 1) {
        if (verbose) {
          message(paste("skip single exon gene", gid))
        }
        next
      }
      if (mat$niso < niso_min) {
        if (verbose) {
          message(paste("skip gene", gid, "for number of iso less than", niso_min))
        }
        next
      }

      #control_ok = FALSE
      #case_ok = FALSE
      #mat$fpkm_mat, rows are samples and columns are isoforms
      # for(k in 1:col(mat$fpkm_mat)) {
      #   if (rownames(mat$fpkm_mat)[k] %in% cases) {
      #     if (all(mat$fpkm_mat[k, ] > FPKM_filter)) {
      #       case_ok = TRUE
      #     }
      #   } else {
      #     if (all(mat$fpkm_mat[k,] > FPKM_filter)) {
      #       control_ok = TRUE
      #     }
      #   }
      # }
      # if (!control_ok || !case_ok) {
      #   if (verbose) {
      #     message(paste("skip gene", gid, "by FPKM filter at", FPKM_filter))
      #   }
      #   next
      # }

      is_ok = FALSE
      for(k in 1:col(mat$fpkm_mat)) {
        if (all(mat$fpkm_mat[, k] > FPKM_filter)) {
          is_ok = TRUE
        }
      }
      if (!is_ok) {
        if (verbose) {
          message(paste("skip gene", gid, "by FPKM filter at", FPKM_filter))
        }
        next
      }

      initpi$g2t_map[[gid]] %in% colnames(initpi$pi)
      gene_pi = initpi$pi[,initpi$g2t_map[[gid]]]
      mlr.data = GenewiseHiddenTable(mat, pi = gene_pi, biasfit=bias_correcter_, normalize = NORMALIZE)
      #mlr.data
      if (is.null(mlr.data)) {
        if (verbose) {
          message(paste("cannot get hidden table. skip gene", gid))
        }
        next
      }
      if (is.null(Y)) {
        gene_ids = rep(gid, ncol(mlr.data$y_mul))
        Y = mlr.data$y_mul
        obs_table = mlr.data$obs_table
        pi = gene_pi
        theta = mlr.data$theta
      } else {
        gene_ids = c(gene_ids, rep(gid, ncol(mlr.data$y_mul)))
        Y = cbind(Y, mlr.data$y_mul)
        obs_table = cbind(obs_table, mlr.data$obs_table)
        pi = cbind(pi, gene_pi)
        new_theta = mlr.data$theta
        lower_half_theta = cbind(matrix(0, nrow = nrow(new_theta), ncol = ncol(theta)), new_theta)
        dim(lower_half_theta)
        upper_half_theta = cbind(theta, matrix(0, nrow = nrow(theta), ncol = ncol(new_theta)))
        theta = rbind(upper_half_theta, lower_half_theta)
      }
      if (is.null(X)) {
        X = mlr.data$X_mul
      }
    }# end if skip_gene
  }
  if (is.null(X)) {
    return (NULL)
  }
  designX = cbind(rep(1,nrow(X)),X)
  colnames(designX)[1] = "const"
  list(pi = pi, gene_ids = gene_ids, X = designX, Y=Y, obs_table = obs_table, theta = theta)
}


#' @importFrom rstan cpp_object_initializer
driver <- function(gene_rf, cases, initpi, priors, bias_coef = bias_coef, BAYES = TRUE, MCMC = FALSE, MAX_ITER = 50, target_genes = NULL, FPKM_filter = 1.0, NORMALIZE = F, verbose=FALSE) {
  result = NULL
  opt_par = NULL
  new_pi = NULL
  new_theta = NULL
  FIX_THETA = TRUE
  if (is.null(bias_coef)) {
    FIX_THETA = FALSE
  }
  result = prepareData(gene_rf, cases, initpi, BAYES, bias_coef, NORMALIZE = NORMALIZE, target_genes=target_genes, FPKM_filter = FPKM_filter, verbose = verbose)
  if (is.null(result)) {
    return (NULL)
  }
  niso = ncol(result$Y)
  ndata = nrow(result$Y)
  Y = result$Y
  X = result$X
  stopifnot(ndata == nrow(X) && ndata == ncol(X))
  OBS = result$obs_table
  theta = result$theta
  pi = result$pi
  fits = NULL
  for (i in 1:MAX_ITER) {
    if (BAYES) {
      HTable = getHiddenTable(OBS, pi, theta)
      dat = list(J = ndata, K = niso, Y = getY(HTable), X = X, a0 = priors$a0, b0 = priors$b0, a1 = priors$a1, b1 = priors$b1, a2 = priors$a2, b2 = priors$b2)
      ##stan est
      #fit = stan(model_code = stan_model_str, model_name = "lca", data = dat, iter = 2012, chains = 4, verbose = TRUE)
      #plot(fit, pars = paste0("beta[",1:ncol(dat$Y),",2]"))
      #traceplot(fit, pars = c("beta"))
      ##stan optimization
      #opt_par = optimizing(stanm, data = dat, hessian = TRUE, algorithm="Newton")

      if (MCMC) {
        invisible(capture.output(stan_result <- rstan::stan(model_code = stan_model_str, data=dat, iter=1000, chains=4, refresh=0)))
        mcmc_pi_mat = as.matrix(summary(stan_result, pars= c("pi"), probs=c(0.5))$summary)
        new_pi = matrix(mcmc_pi_mat[,1], byrow=TRUE, nrow=nrow(result$Y), ncol=ncol(result$Y), dimnames = list(rownames(result$Y), colnames(result$Y)))
      } else {
        invisible(capture.output(opt_par <- rstan::optimizing(stanm, data = dat, hessian = TRUE, algorithm="LBFGS", seed=17)))
        new_pi = matrix(opt_par$par[grepl("pi", names(opt_par$par))], nrow=nrow(result$Y), ncol=ncol(result$Y), dimnames = list(rownames(result$Y), colnames(result$Y)))
      }

      if (!FIX_THETA) {
        new_theta = getTheta(getHiddenTable(OBS, new_pi, theta))
      }
    } else {
      HTable = getHiddenTable(OBS, pi, theta)
      Y = getY(HTable)
      fits <- nnet::multinom(getY(HTable) ~ X - 1, trace = FALSE, MaxNWts = 100000)
      new_pi = fitted(fits)
      rownames(new_pi) = rownames(result$Y)
      colnames(new_pi) = colnames(result$Y)

      if (!FIX_THETA) {
        new_theta = getTheta(getHiddenTable(OBS, new_pi, theta))
      }
    }
    if ( !is.null(new_pi)) {
      if(norm(new_pi - pi) < 0.001) {
        if(verbose) {
          print(paste("terminated at iter", i))
        }
        break
      }
    }
    pi = new_pi
    if (!FIX_THETA) {
      theta = new_theta
    }
  }

  pvalues = NULL
  betas = NULL
  if (BAYES) {
    if (MCMC){
      beta_mat = as.matrix(summary(stan_result, pars=c("beta"))$summary)
      beta_idx = grepl("beta\\[.*,2\\]", rownames(beta_mat))
      betas = beta_mat[beta_idx,1]
      z = beta_mat[beta_idx,1] / beta_mat[beta_idx, 2]
      pvalues <- (1 - pnorm(abs(z), 0, 1)) * 2
    } else {
      beta_idx = grepl("beta\\[.*,2\\]", names(opt_par$par))
      betas = opt_par$par[beta_idx]
      #inv_h = solve(opt_par$hessian[which(beta_idx), which(beta_idx)])
      #se = sqrt(-diag(inv_h))
      pvalues = tryCatch({
        se = sqrt(diag(solve(-opt_par$hessian)[which(beta_idx),which(beta_idx)]))
        if (sum(is.na(se))) {
          invisible(capture.output(stan_result <- stan(model_code = stan_model_str, data=dat, iter=1000, chains=4, refresh=0, seed=42)))
          beta_mat = as.matrix(summary(stan_result, pars=c("beta"))$summary)
          beta_idx = grepl("beta\\[.*,2\\]", rownames(beta_mat))
          betas = beta_mat[beta_idx,1]
          z = beta_mat[beta_idx,1] / beta_mat[beta_idx, 2]
          (1 - pnorm(abs(z), 0, 1)) * 2
        } else {
          z = betas/se
          (1 - pnorm(abs(z), 0, 1)) * 2
        }
      }, error = function(e) {
        message("hessian not inversible, trying MCMC")
        #invisible(capture.output(stan_result <- stan(model_code = stan_model_str, data=dat, iter=1000, chains=4, refresh=0, seed = 42)))
        #beta_mat = as.matrix(summary(stan_result, pars=c("beta"))$summary)
        #beta_idx = grepl("beta\\[.*,2\\]", rownames(beta_mat))
        #betas = beta_mat[beta_idx,1]
        #z = beta_mat[beta_idx,1] / beta_mat[beta_idx, 2]
        #(1 - pnorm(abs(z), 0, 1)) * 2
        rep(1.0, length(betas))
      })
      if (verbose) {
        print(paste("pvalue", pvalues))
      }
    }
  } else {
    #print(summary(fits))
    betas = summary(fits)$coefficients
    z <- summary(fits)$coefficients/summary(fits)$standard.errors
    pvalues <- (1 - pnorm(abs(z), 0, 1)) * 2
  }
  list(pvalues = pvalues, gene_ids = result$gene_ids, pi = pi, betas = betas, niso=niso)
}

####
####
####



main <-function(gene_rf, cases, BAYES, priors, MAX_ITER, genelist = NULL, start = 1, FPKM_filter = 1.0, STEP = 1, RANDOM_PI_START = F, BIAS=F, NORMALIZE = F, MCMC=FALSE, verbose=F){
  if (is.null(genelist)) {
    genelist = names(gene_rf[[1]])
  }
  nblock = ceiling(length(genelist) / STEP)
  pb <- txtProgressBar(min = 0, max = nblock, style = 3)

  bias_coef = NULL
  if (BIAS) {
    print("using bias estimation")
    total_bias_mat = preprocess(gene_rf)
    bias_coef = getBiasCoef(total_bias_mat)
  }
  pi = GetPiMatrix(gene_rf, RANDOM_PI_START)
  gene_rf[[1]][[1]]
  collector = list()
  for (i in 1:nblock) {
    # update progress bar
    setTxtProgressBar(pb, i)
    e = start + STEP - 1
    targets = SubsetGeneIds(genelist, s=start, e=e)
    if (is.null(targets)) next
    if (verbose) {
      print("==========================")
      print(targets)
    }
    result = driver(gene_rf, cases, pi, bias_coef, BAYES=BAYES, priors = priors, FPKM_filter = FPKM_filter, target_genes = targets, MAX_ITER=MAX_ITER, MCMC = MCMC, NORMALIZE = NORMALIZE, verbose = verbose)
    if (!is.null(result)) {
      collector[[length(collector)+1]] = result
      names(collector)[length(collector)] = unique(result$gene_ids)
    }
    start = e + 1
  }
  collector
}

#' Load data from strawberry input.
#'
#' @param dir A directory where the strawberry exonpath context output are stored.
#' @return A list of of samples in Raspberry data format.
#' @examples loadData("inst/extdata/AT_RD100_exonpath_context/")
#' @export
loadData <- function(dir) {
  files = list.files(path = dir)
  file_paths = paste0(dir, "/", files)
  nfiles = length(files)
  gene_rf = vector("list", length = nfiles)
  for (i in 1:nfiles) {
    print(paste("loading file", file_paths[i]))
    rf =  read.table(file_paths[i], header=T, stringsAsFactors = F, sep = "\t")
    split_rf = split(rf, rf$gene_id)
    gene_rf[[i]] = split_rf
  }
  gene_rf
}

#' Run differential splicing anaysis.
#'
#' @param raspberry_df A list of samples in Raspberry input data format
#' @param cases A vector of strings indicating which samples are treatment samples.
#' @param MAX_ITER A number for the max iteration num for the EM algorithm. Default is 8.
#' @param BIAS A boolean for doing bias correction or not, default is False
#' @param FPKM_filter A number for fpkm filtering. Default is 1.0.
#' @param MCMC A boolean for using MCMC instead of EM algorithm for optimization. Note that this will be very slow. Default is False.
#' @param verbose A boolean for more verbose messages. Default is False.
#' @return A list of genes. Each gene has multiple data.frame
#' @examples
#' diffiso(raspdata, c("caseSample1, caseSample2, caseSample3"), MAX_ITER = 5, BIAS = F, FPKM_filter = 1.0, MCMC = FALSE, verbose = FALSE)
#' @export
diffiso <- function(raspberry_df, cases, MAX_ITER = 8, BIAS = F, FPKM_filter = 0.1, MCMC = FALSE, verbose = FALSE) {
  print("Estimating empirical bayes priors")
  collector1 = main(raspberry_df, cases, BAYES = FALSE, MAX_ITER = 3, BIAS = BIAS, FPKM_filter = FPKM_filter, RANDOM_PI_START = T, verbose = verbose)

  priors = getPrior(collector1, nsamples = nrow(collector1[[1]]$pi), trim=0.25)

  print("Calculating differential alternative spliced transcripts")
  collector2 =  main(raspberry_df, cases, BAYES = TRUE, priors = priors, MAX_ITER = MAX_ITER, BIAS = BIAS, MCMC = MCMC, FPKM_filter = FPKM_filter,
                     RANDOM_PI_START = T)
  collector2
}


#' Run differential splicing anaysis.
#'
#' @param dataDir A string for the path of the input data dir.
#' @param cases A vector of strings indicating which samples are treatment samples.
#' @param MAX_ITER A number for the max iteration num for the EM algorithm. Default is 8.
#' @param BIAS A boolean for doing bias correction or not, default is False
#' @param FPKM_filter A number for fpkm filtering. Default is 1.0.
#' @param MCMC A boolean for using MCMC instead of EM algorithm for optimization. Note that this will be very slow. Default is False.
#' @param verbose A boolean for more verbose messages. Default is False.
#' @return A list of genes. Each gene has multiple data.frame
#' @examples
#' run_analysis("my_path_to_data", c("caseSample1, caseSample2, caseSample3"), MAX_ITER = 5, BIAS = F, FPKM_filter = 1.0, MCMC = FALSE, verbose = FALSE)
#' @export
run_analysis <- function(dataDir, cases, MAX_ITER = 8, BIAS = F, FPKM_filter = 0.1, MCMC = FALSE, verbose = FALSE) {
  gene_rf = loadData(dir = dataDir)
  list(result=diffiso(gene_rf, cases, MAX_ITER, BIAS, FPKM_filter, MCMC, verbose))
}
