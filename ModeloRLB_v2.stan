data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  int<lower=0,upper=1> y[N];
  vector<lower=0>[P] alpha_pi;
  vector<lower=0>[P] beta_pi;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_alpha;
}

parameters {
  real B_0;
  vector[P] B_j;
  vector<lower=0,upper=1>[P] pi_j;
}

transformed parameters {
  // Coeficientes 
  vector[P] effective_beta = B_j .* pi_j;
}

model {
  // Priors 
  B_0 ~ normal(0, sigma_alpha);
  B_j ~ normal(0, sigma_beta);
  pi_j ~ beta(alpha_pi, beta_pi);
  
  // Likelihood 
  y ~ bernoulli_logit(B_0 + X * effective_beta);
}

generated quantities {
  // Solo lo ESENCIAL para model comparison
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | B_0 + dot_product(X[n], effective_beta));
  }
}

