data {
  int<lower=1>   T;             // years
  vector[T]      flow;          // mean flow each trapping season
  int<lower=0>   released[T];   // # fish released in efficiency trials
  int<lower=0>   caught[T];     // # recaptured in trials
  int<lower=0>   C[T];          // screw‐trap catch each year
  int<lower=0>   E_obs[T];      // escapement counts each year
}
parameters {
  // -- trap detection
  real          gamma0;         // intercept
  real          gamma1;         // slope on flow

  // -- stock–recruit
  real<lower=0> a;              // BH productivity
  real<lower=0> K;              // habitat capacity
  real<lower=0> sigma_sr;       // observation / process error

  // latent true fry for each year
  vector<lower=0>[T] N_fry;
}
model {
  // 1) Priors
  gamma0   ~ normal(0, 5);
  gamma1   ~ normal(0, 1);
  a        ~ normal(0.5, 1);
  K        ~ normal(50000, 30000);
  sigma_sr ~ normal(0, 5000);

  // 2) Trap‐efficiency trials → identify gamma0,gamma1
  for (t in 1:T) {
    caught[t] ~ binomial(released[t],
                         inv_logit(gamma0 + gamma1 * flow[t]));
  }

  // 3) True fry → link to catch
  for (t in 1:T) {
    real p_tr = inv_logit(gamma0 + gamma1 * flow[t]);
    // N_fry is latent; catch ~ Binomial(N_fry, p_tr)
    C[t]    ~ binomial(N_fry[t], p_tr);
  }

  // 4) Stock–recruit (Beverton–Holt) → escapement
  for (t in 1:T) {
    real muE = a * N_fry[t] / (1 + N_fry[t] / K);
    E_obs[t] ~ normal(muE, sigma_sr);
  }
}
