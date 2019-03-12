mytdist <- function(x, m, s, df) {
  lambda = 1e-6
  regp = log(dt( (x-m)/s, df)/s) - lambda * sum(m^2 + s^2 + df^2)
  exp(regp)
}

getPrior <- function(collector, nsamples, trim = 0.1) {

  betas = NULL
  for(i in 1:length(collector)) {
    if (i == 1) {
      betas = collector[[i]]$betas
    }
    else {
      betas = rbind(betas, collector[[i]]$betas)
    }
  }
  a = c()
  b = c()
  estimates = vector("list", length = nsamples)
  for(j in 1:nsamples) {
     tmp = betas[, j][betas[,j] > quantile(betas[,j], trim) & betas[,j] < quantile(betas[,j], 1-trim) ]
     #hist(a, breaks=30)
     res= MASS::fitdistr(tmp, mytdist, list(m=0, s = 0.02, df = 20))
     a_j = res$estimate[3]/2
     b_j = (res$estimate[2] * res$estimate[2] * res$estimate[3]) / 2
     a_j
     b_j
     a = c(a, a_j)
     b = c(b, b_j)
     estimates[[j]] = res$estimate
  }
  list(estimates = estimates, a = a, b = b, betas = betas, a0 = a[1], a1 = a[2], a2 = a[3], b0 = b[1], b1 = b[2], b2 = b[3])
  #list(a0 = alpha[1], a1 = mean(alpha[2:5]), b0 = beta[1], b1 = mean(beta[2:5]))
}

winsor.var<-function(x, trim=0.1) {
  xq<-quantile(x,probs=c(trim, 1 - trim))
  x[x < xq[1]]<-xq[1]
  x[x > xq[2]]<-xq[2]
  return(var(x))
}

winsor.sd <- function(x, trim=0.1) {
  sqrt(winsor.var(x, trim))
}

winsorize <-function(x, trim=0.1) {
  xq<-quantile(x,probs=c(trim, 1 - trim))
  x[x < xq[1]]<-xq[1]
  x[x > xq[2]]<-xq[2]
  return(x)
}
#getPrior <- function(collector, lq = 0.025, uq = 0.975) {
#  taus = getPrecision(collector)
#  quant=quantile(taus, c(lq, uq))
#  fitgamma = fitdistr(taus[taus<quant[2] & taus > quant[1]], dgamma, start=list(shape=25, rate=50))
#  fitgamma$estimate
#}
?fitdistr
