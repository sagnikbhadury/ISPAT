
#----------------------------
# extra function
#----------------------------

ext_loop <- function(count, c, N_c, pos, Y, S, Clusters, KS_c, lS, stan_model_expo, stan_model_matern, sGLM_method, Kernel){
  storeZ <- NA
  mydata <- Y
  pos <- which(S[,3] == Clusters[c]); S_loc_c = S[pos,-3]
  stan_Data = list(N=N_c[c],y=as.vector(as.matrix(mydata)[count,pos]),Sx=S_loc_c[,1],Sy=S_loc_c[,2],lS=lS) ## Input data for stan models
  if(Kernel == "RBF"){ ## For Kernel choice of exponential (RBF) kernel
    if(sGLM_method == "MLE"){ ## For MLE estimation of sGLM method with RBF Kernel
      MLE_RBF_G <- rstan::sampling(stan_model_expo, data = stan_Data, iter = 500, chains = 2, warmup = 200, thin = 1,verbose = TRUE)
      MLEs <- summary(MLE_RBF_G)$summary[,"mean"]
      beta_est <- MLEs[1] ## MLE of Beta
      sigmaS_est <- MLEs[2] ## MLE of sigmaS
      sigmaEPS_est <- MLEs[3]  ## MLE of sigmaEPS
    } ## if loop end for sGLM_method == "MLE" with RBF Kernel
    if(sGLM_method == "VB"){
      ## For VB estimation of sGLM method with RBF Kernel
      print("We recommend 'sampling' from the posterior for sGLM!")
      VB_RBF_G <- rstan::vb(stan_model_expo, data = stan_Data, seed = 50, output_samples = 5000, refresh = 0)
      VBests <- summary(VB_RBF_G)$summary[,"mean"]
      beta_est <- VBests[1] ## VBE of Beta
      sigmaS_est <- VBests[2] ## VBE of sigmaS
      sigmaEPS_est <- VBests[3] ## VBE of sigmaEPS
    } ## if loop end for sGLM_method == "VB" with RBF Kernel
  } ## if loop end for Kernel = "RBF"
  if(Kernel == "Matern"){ ## For Kernel choice of Matern kernel
    if(sGLM_method == "MLE"){ ## For MLE estimation of sGLM method with Matern Kernel
      MLE_Matern_G <- rstan::sampling(stan_model_matern, data = stan_Data, iter = 500, chains = 2, warmup = 200, thin = 1,verbose = TRUE)
      MLEs <- summary(MLE_Matern_G)$summary[,"mean"]
      beta_est <- MLEs[1] ## MLE of Beta
      sigmaS_est <- MLEs[2] ## MLE of sigmaS
      sigmaEPS_est <- MLEs[3]  ## MLE of sigmaEPS
    } ## if loop end for sGLM_method == "MLE" with Matern Kernel
    if(sGLM_method == "VB"){ ## For VB estimation of sGLM method with Matern Kernel
      # Run variational inference
      print("We recommend 'sampling' from the posterior for sGLM!")
      VB_Matern_G <- rstan::vb(stan_model_matern, data = stan_Data, seed = 50, output_samples = 5000, refresh = 0)
      VBests <- summary(VB_Matern_G)$summary[,"mean"]
      beta_est <- VBests[1] ## VBE of Beta
      sigmaS_est <- VBests[2] ## VBE of sigmaS
      sigmaEPS_est <- VBests[3] ## VBE of sigmaEPS
    } ## if loop end for sGLM_method == "VB" with RBF Kernel
  } ## if loop end for Kernel = "Matern"

  if(sigmaS_est< 0.01 || sigmaEPS_est< 0.01){
    V <- (sigmaS_est+0.001)*KS_c[[c]] + (sigmaEPS_est+0.001)*diag(N_c[c])
  }else{
    V <- (sigmaS_est)*KS_c[[c]] + (sigmaEPS_est)*diag(N_c[c])
  }
  #####   Calculation of the latent gene expression matrix    #####
  if(sigmaEPS_est < 0.001){
    storeZ <- ((sigmaEPS_est+0.001)*solve(V))%*%(as.matrix(Y)[count,pos]-beta_est)
  }else{
    storeZ <- ((sigmaEPS_est)*solve(V))%*%(as.matrix(Y)[count,pos]-beta_est)
  }
  rm(mydata)
  res<- data.frame(storeZ)
  return(res)
} ## End of for each gene loop
