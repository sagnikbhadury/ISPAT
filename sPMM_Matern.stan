// Spatial Poisson mixed model (sPMM) with Matern kernel

data {
  int<lower=0> N;             // Number of data points
  array[N] int<lower=0> counts; //int<lower=0> counts[N];     //array[N] int<lower=0> counts;    // Observed counts
  array[N] real Sx;   //real Sx[N];         // Spatial locations (x co-ordinates)
  array[N] real Sy;   //real Sy[N];         // Spatial locations (y co-ordinates)
  real<lower=0> lS;   //length scale parameter of the spatial kernel (KS)
}

parameters {
  real beta;               // Fixed effects coefficients
  vector[N] eta;                // Spatial random effects
  real<lower=0> sigmaS;   //scale parameter of GP for time function f_{S}()
  real<lower=0> sigmaEPS; //Error variance
}

transformed parameters {
  vector[N] lambda = exp(beta + eta); // Mean of the Poisson (log link)
}

model {
  
  matrix[N, N] KS = (gp_matern32_cov(Sx,sigmaS,lS).*gp_matern32_cov(Sy,sigmaS,lS))/(sigmaS^2);

    // diagonal elements
    for (n in 1:N) {
      KS[n, n] = KS[n, n] + square(sigmaEPS);
    }

  eta ~ multi_normal(rep_vector(0, N), KS);
  
  // Priors
  beta ~ std_normal();
  sigmaS ~ std_normal(); 
  sigmaEPS ~ std_normal();
  
  // Likelihood
  counts ~ poisson_log(lambda);
}
