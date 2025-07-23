#------------------------------------------------------------------------------------
# Gaussian Data input for VB ISPAT model with clusters and reference prior
#------------------------------------------------------------------------------------


pacman::p_load(doParallel, foreach, rstan, Matrix, matrixcalc, BiocManager, cmdstanr, parallel, CVglasso)
#----------------------

source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/domain informed/svi_msfa.R")
source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/domain informed/cavi_msfa.R")
source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/domain informed/stan programs.R")

#-----------------------
# specify data type
#-----------------------

data_type = "Gaussian"

if(data_type == "Gaussian"){
  stan_model_matern <- stan_model(model_code = mystan_gauss_matern)
  stan_model_expo <- stan_model(model_code = mystan_gauss_expo)
}



#----------------------------------------------
# gaussian SpaceX model with Reference Prior
#----------------------------------------------


ISPAT <- function(Y,S,ncores, RefPrior,
                                   use_ref = c(TRUE,FALSE),
                                   Kernel = c("RBF","Matern"),
                                   sGLM_method = c("MLE","VB"),
                                   VB_MSFA = c(TRUE, FALSE),
                                   MSFA_method = c("CAVI","SVI")){

  ## Inputs

  ## Y: GxN Gene expression matrix (N: no of spatial locations, G: no of genes)
  ## S: Spatial location matrix (Nx3 matrix), First two columns consists of spatial locatiosn and last column has cluster annotations
  ## ncores: number of cores to run in parallel
  ## use_ref: Use reference prior for the Multi-study Factor Analysis TRUE/FALSE
  ## RefPrior: Reference prior for the Multi-study Factor Analysis : USE SQUARED of ACTUAL COV MATRIX
  ## Kernel: choice of kernels for the spatial parameter are RBF and Matern
  ## sGLM_method: Estimation method based on MLE or VB
  ## MSFA_method: MSFA method can be implemented with two different type of algorithms (CAVI ro SVI)

  ## SVI not giving good results, use CAVI
  # generated_data <- readRDS(paste0("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/domain informed/simulations/maySim/simData","_N_1_",N_c[1], "_N_2_",N_c[2] , "_N_3_", N_c[3],"_cells_",G,"_repli_",5,"_",Matern,"_ls50_datagenerate.rds"), refhook = NULL)
  # Y <- generated_data$Obs_matrix[[1]]
  # S <- generated_data$Coords_annotations
  # RefPrior = diag(1, nrow = dim(Y)[1], ncol = dim(Y)[1])
  # ncores = detectCores() -3
  # use_ref = FALSE ; Kernel = "Matern" ; sGLM_method = "VB"; VB_MSFA = TRUE; MSFA_method = "CAVI"

  ######   Global Parameters   ########
  N <- dim(Y)[2]  # Number of spatial locations
  G <- dim(Y)[1]  # Number of genes
  C <- length(unique(S[,3])) # Number of clusters
  Clusters <- unique(S[,3])

  N_c <- numeric()
  sigmaS_est <- matrix(0,G,C)
  sigmaEPS_est <- matrix(0,G,C)
  Z_est <- list() ##latent gene expression matrix
  KS_c <- list()  ##Cluster specific spatial kernels

  ### Cluster sizes ####
  for (c in 1:C) {
    pos <- which(S[,3] == Clusters[c])
    N_c[c] <- length(pos)
    Z_est[[c]] <- matrix(0,nrow = N_c[c],ncol = G)
    KS_c[[c]] <- matrix(0,N_c[c],N_c[c])
  }

  for (c in 1:C){
    pos <- which(S[,3] == Clusters[c])
    ### Rho estimation ###
    a <- dist(S[pos,-3])
    a_max <- log10(2*max(a))
    a_min <- log10(min(a)/2)
    a_seq <- seq(a_min, a_max, length.out = 10)
    lS <- 10^(a_seq[5])

    if(Kernel == "RBF"){
      S_loc_c = S[pos,-3]
      for (i in 1:N_c[c]) {
        for (j in 1:i) {
          dist_loc <- (S_loc_c[i,1] - S_loc_c[j,1])^2 + (S_loc_c[i,2] - S_loc_c[j,2])^2
          KS_c[[c]][i,j] <- exp(-dist_loc/(2*lS^2))
          KS_c[[c]][j,i] <- KS_c[[c]][i,j] }}} #for symmetry

    if(Kernel == "Matern"){
      S_loc_c = S[pos,-3]
      for (i in 1:N_c[c]) {
        for (j in 1:i) {
          distS <- sqrt((S_loc_c[i,1] - S_loc_c[j,1])^2 + (S_loc_c[i,2] - S_loc_c[j,2])^2)/lS
          KS_c[[c]][i,j] <- (1+distS*sqrt(3))*exp(-distS*sqrt(3))
          KS_c[[c]][j,i] <- KS_c[[c]][i,j] }}
      if(is.positive.definite(KS_c[[c]])==FALSE){ KS_c[[c]] <- as.matrix(nearPD(KS_c[[c]])$mat)  }}

    print(paste0("Fitting Parallel Spatial Mixed Effects Model for cluster - ", c))
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    Z_est[[c]]<- suppressWarnings({ foreach(my_count = 1:G, .packages = c('base','Matrix', 'rstan', 'matrixcalc','CVglasso'), .combine = cbind,
                .export = ls(globalenv())) %dopar% {
                 myresult <- ext_loop(my_count, c, N_c, pos, Y, S, Clusters, KS_c, lS, stan_model_expo, stan_model_matern, sGLM_method, Kernel) }})
    ## Stop the cluster
    stopCluster(cl)
    print(paste0("Finished with Parallel Spatial Mixed Effects Model for cluster - ", c))
  }
  print("Finished with Parallel Spatial Mixed Effects Model for all clusters")

  refs <- lapply(1:C, function(c) { Glasso <- CVglasso::CVglasso( X = Z_est[[c]], nlam = 20, lam.min.ratio = 0.0000001, diagonal = FALSE,
                      crit.cv = "BIC", maxit = 120, adjmaxit = NULL, K = 10, cores = 10); Glasso <- Glasso$Omega ; #Glasso <- abs(Glasso); Glasso[Glasso < 0.01] <- 0; Glasso <- as.matrix(Glasso)
                      }); #for (f in 1:length(N_c)) {diag(refs[[f]])<-0};

  ####   Multi-study Factor model on latent gene expression matrix   ####

  print("fitting Multi-Study Factor Model")
  if (VB_MSFA == TRUE) {
    if (MSFA_method == "CAVI") {
      VBfit_MSFA <- cavi_msfa(Z_est, floor(2 * log(G)), rep(floor(2 * log(G)), C), scale =  FALSE)
      Shared_Net <- tcrossprod(VBfit_MSFA$mean_phi)
      Cluster_Net <-  lapply(1:C, function(c) { tcrossprod(VBfit_MSFA$mean_lambda_s[[c]]) + tcrossprod(VBfit_MSFA$mean_phi) })
      # ref_net <- Cluster_Net
      if(use_ref == TRUE){
        ref_net <- lapply(1:C, function(c) { 1/sqrt(refs[[c]]) %*% Cluster_Net[[c]] %*% 1/sqrt(refs[[c]]) })
      }
    } else if (MSFA_method == "SVI") {
      VBfit_MSFA <- svi_msfa(Z_est, floor(2 * log(G)), rep(floor(2 * log(G)), C), scale = FALSE)
      Shared_Net <- tcrossprod(VBfit_MSFA$mean_phi)
      Cluster_Net <-  lapply(1:C, function(c) { tcrossprod(VBfit_MSFA$mean_lambda_s[[c]]) }) }
      # ref_net <- Cluster_Net
      if(use_ref == TRUE){
      ref_net <- lapply(1:C, function(c) { (1/sqrt(refs[[c]])) %*% Cluster_Net[[c]] %*% (1/sqrt(refs[[c]])) })
      }
    } else {
    print("Please set VB_MSFA = TRUE and rerun"); break;
    }
  print("Function Finished - Happy Exploring!")
  return(list(Shared_Net=Shared_Net,Cluster_Net=Cluster_Net, glasso_refs = refs ,KS_c = KS_c, VBfit_MSFA = VBfit_MSFA, Z_est = Z_est))
}

