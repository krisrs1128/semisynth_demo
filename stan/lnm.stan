data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> D;
  matrix[N, D] x;
  array[N, K + 1] int<lower=0> y;
}

parameters {
  matrix[D, K] B;
  vector<lower=0>[K] sigmas;
  matrix[N, K] mu;
  real<lower=0> sigma2_beta;
}

model {
  vector[K + 1] p;
  vector[K] e_mu;
  for (i in 1:N) {
    mu[i] ~ multi_normal(B' * to_vector(x[i]), diag_matrix(sigmas));
    e_mu = to_vector(exp(mu[i]));
    p = append_row(e_mu / (1 + sum(e_mu)), 1 / (1 + sum(e_mu)));
    y[i] ~ multinomial(to_vector(p));
  }
}
