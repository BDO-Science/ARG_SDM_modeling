data {
  int<lower=1> J;              
  int<lower=1> Kb;             
  matrix[J,Kb] B;              
  vector[J] X;                 
  real<lower=0> logN_max;      
  int<lower=0> n1[J];          // marked released
  int<lower=0> m2[J];          // marked recaptured
  int<lower=0> u2[J];          // unmarked catch
  // optionally: vector[J] flow_effort; // if you want p covariates
}
parameters {
  vector[Kb] gamma;            
  real       phi;              
  vector[J]  v;                
  real<lower=0> sigma_v;       
  vector[J]  logU;             
  vector[1]  beta_p;           // if you want constant p or covariate on p
}
model {
  // –– abundance spline + covariate + v
  for (k in 3:Kb)
    (gamma[k] - 2*gamma[k-1] + gamma[k-2]) ~ normal(0, 1);
  gamma[1:2] ~ normal(0,10);
  v          ~ normal(0, sigma_v);
  sigma_v    ~ gamma(1,0.05);
  phi        ~ normal(0,5);

  for (t in 1:J) {
    // build log‐U with truncation
    real lu = dot_product(B[t], gamma) + phi * X[t] + v[t];
    lu      = fmin(lu, logN_max);  
    // capture prob (simple constant here; or add covariates)
    real p = inv_logit(beta_p[1]);
    //  exact binomial on marks:
    m2[t] ~ binomial(n1[t], p);
    //  *if* we accept the Poisson approx for the unmarked:
    u2[t] ~ poisson(exp(lu) * p);
    // else you’d need discrete U and
    //   u2[t] ~ binomial(round(exp(lu)), p);
  }
}
generated quantities {
  real Utot = sum(exp(logU));  
  real p_est = inv_logit(beta_p[1]);
}
