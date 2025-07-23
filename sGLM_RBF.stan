// Spatial Gaussian linear mixed model (sGLM) with exponential kernel

data {
  int<lower=0> N;     // No of spatial locations (N) 
  vector[N] y;        // Normalized Gene expression vector (N x 1) for a particular gene
  array[N] real Sx;         //real Sx[N];          // Spatial locations (x co-ordinates)
  array[N] real Sy;         //real Sy[N];          // Spatial locations (y co-ordinates)
  real<lower=0> lS;   //length scale parameter of the spatial kernel (KS)
}

parameters {
  real beta;              //intercept
  real<lower=0> sigmaS;   //scale parameter of GP for time function f_{S}()
  real<lower=0> sigmaEPS; //Error variance
}


model {
   matrix[N, N] L_KS;
   
   matrix[N, N] KS = (gp_exp_quad_cov(Sx,sigmaS,lS).*gp_exp_quad_cov(Sy,sigmaS,lS))/(sigmaS^2);

    // diagonal elements
    for (n in 1:N) {
      KS[n, n] = KS[n, n] + square(sigmaEPS);
    }
  
   L_KS = cholesky_decompose(KS);
   
   //Priors
   beta ~ std_normal();
   sigmaS ~ std_normal(); 
   sigmaEPS ~ std_normal(); 
   
  //Model
  y ~ multi_normal_cholesky(rep_vector(beta, N), L_KS);
}


