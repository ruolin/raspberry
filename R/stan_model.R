library(rstan)
stan_model_str <- "
  data {
    int<lower = 2> J; // num samples
    int<lower = 1> K; // num transcripts
    vector<lower = 0>[K] Y[J]; // counts
    vector[J] X[J]; //predictors
    real<lower = 0> a0; // fixed prior of const effect, shape
    real<lower = 0> b0; // fixed prior of const effect, rate
    real<lower = 0> a1; // fixed prior of condition effect, shape
    real<lower = 0> b1; // fixed prior of condition effect, rate
    real<lower = 0> a2; // fixed prior of sample effect, shape
    real<lower = 0> b2; // fixed prior of sample effect, rate
  }

  parameters {
    matrix[K,J] beta;
    vector<lower = 0> [J] gamma;
    #real<lower = 0> gamma_c;
  }

  transformed parameters {
    simplex[K] pi[J];
    vector<lower = 0> [J] sigma;

    for (j in 1:J) {
      pi[j] = softmax(beta * X[j]);
    }

    #sigma[1] = 1 / sqrt(gamma_c);
    #sigma[1] = 1/gamma_c;

    for (j in 1:(J)) {
      sigma[j] = 1 / sqrt(gamma[j]);
      #sigma[j] = 1 / gamma[j];
    }
  }

  model {
    //priors
    for (j in 1:J) {
      if (j == 1) {
        gamma[j] ~ gamma(a0, b0);
      } else if (j == 2) {
        gamma[j] ~ gamma(a1, b1);
      } else {
        gamma[j] ~ gamma(a2, b2);
      }
    }
    #gamma ~ gamma(a0, b0);
    for (j in 1:(J)) {
      col(beta,j) ~ normal(0, sigma[j]);
    }

    for (j in 1:J) {
      for (k in 1:K) {
        target += Y[j][k] * log (pi[j][k]);
      }
    }
  }
"


options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)
stanm = rstan::stan_model(model_code = stan_model_str, model_name= model_name)

